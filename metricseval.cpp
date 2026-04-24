#include "metricseval.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <vector>

namespace
{
constexpr double EPS_POWER = 1e-30;
constexpr int PSD_WIN_SIZE = 2048;
constexpr int PSD_OVERLAP  = 1024;
constexpr int PSD_NFFT     = 8192;

static double mW_to_dBm(double p_mW)
{
    if (p_mW <= 0.0)
        return -std::numeric_limits<double>::infinity();

    return 10.0 * std::log10(p_mW);
}

static double dBm_to_mW(double p_dBm)
{
    return std::pow(10.0, p_dBm / 10.0);
}

static double avgPower_mW(const std::vector<std::complex<double>>& x)
{
    if (x.empty())
        return 0.0;

    double sum = 0.0;

    for (const auto& s : x)
        sum += std::norm(s);

    return sum / static_cast<double>(x.size());
}
}

MetricsEval::MetricsEval()
    : fftw8192(FFT_SIZE)
{
    window = hamming(WINDOW_SISE);
}

std::vector<std::complex<double>> MetricsEval::normalizeSignal(
    const std::vector<std::complex<double>>& tx)
{
    if (tx.empty())
        return {};

    double sum = 0.0;

    for (const auto& s : tx)
        sum += std::norm(s);

    double rms = std::sqrt(sum / static_cast<double>(tx.size()));

    if (rms <= 0.0)
        return tx;

    std::vector<std::complex<double>> tx_norm(tx.size());

    for (size_t i = 0; i < tx.size(); ++i)
        tx_norm[i] = tx[i] / rms;

    return tx_norm;
}

double MetricsEval::computeNMSE_dB(const std::vector<std::complex<double>>& ref,
                                   const std::vector<std::complex<double>>& test)
{
    if (ref.empty() || test.empty() || ref.size() != test.size())
        return 0.0;

    std::complex<double> num(0.0, 0.0);
    double den = 0.0;

    for (size_t i = 0; i < ref.size(); ++i)
    {
        num += std::conj(ref[i]) * test[i];
        den += std::norm(ref[i]);
    }

    if (den <= 0.0)
        return 0.0;

    std::complex<double> a = num / den;

    double errPow = 0.0;
    double refPow = 0.0;

    for (size_t i = 0; i < ref.size(); ++i)
    {
        std::complex<double> d = a * ref[i];
        errPow += std::norm(test[i] - d);
        refPow += std::norm(d);
    }

    if (refPow <= 0.0)
        return 0.0;

    return 10.0 * std::log10(errPow / refPow);
}

void MetricsEval::computePSD(
    const std::vector<std::complex<double>>& tx,
    double Fs,
    double oversample,
    std::vector<double>& freq,
    std::vector<double>& psd,
    std::vector<std::complex<double>>& sig_buff)
{
    if (tx.size() < WINDOW_SISE)
    {
        if (sig_buff.size() > 4 * WINDOW_SISE)
            sig_buff.clear();

        for (size_t i = 0; i < tx.size(); ++i)
            sig_buff.push_back(tx[i]);

        if (sig_buff.size() < 4 * WINDOW_SISE)
            return;

        freq.clear();
        psd.clear();
        computePSDWelch(sig_buff, Fs, oversample, freq, psd);
        return;
    }

    freq.clear();
    psd.clear();
    computePSDWelch(tx, Fs, oversample, freq, psd);
}

QPair<double, double> MetricsEval::Calc_BER(std::vector<Symbols>& symbols)
{
    if (symbols.empty() || symbols[0].data_tx.empty())
        return {0.0, 0.0};

    double errors_noDPD = 0.0;
    double errors_withDPD = 0.0;
    double totalBits = 0.0;

    for (int n = 0; n < symbols.size(); ++n)
    {
        const size_t count = std::min({symbols[n].data_tx.size(),
                                       symbols[n].data_rx.size(),
                                       symbols[n].data_rx_with_DPD.size()});

        for (size_t i = 0; i < count; ++i)
        {
            if (symbols[n].data_rx[i] != symbols[n].data_tx[i])
                errors_noDPD++;

            if (symbols[n].data_rx_with_DPD[i] != symbols[n].data_tx[i])
                errors_withDPD++;
        }

        totalBits += static_cast<double>(count);
    }

    if (totalBits <= 0.0)
        return {0.0, 0.0};

    return {errors_noDPD / totalBits, errors_withDPD / totalBits};
}

QPair<double, double> MetricsEval::Calc_EVM(const std::vector<Symbols>& frames)
{
    double err_no_dpd = 0.0;
    double err_with_dpd = 0.0;
    double ref_power = 0.0;

    for (const auto& f : frames)
    {
        size_t n = std::min({f.tr_sym_clean.size(),
                             f.rec_sym_noisy.size(),
                             f.rec_sym_noisy_with_DPD.size()});

        for (size_t i = 0; i < n; ++i)
        {
            const auto& ref = f.tr_sym_clean[i];

            err_no_dpd += std::norm(f.rec_sym_noisy[i] - ref);
            err_with_dpd += std::norm(f.rec_sym_noisy_with_DPD[i] - ref);
            ref_power += std::norm(ref);
        }
    }

    QPair<double, double> res{0.0, 0.0};

    if (ref_power > 0.0)
    {
        res.first = 20.0 * std::log10(std::sqrt(err_no_dpd / ref_power));
        res.second = 20.0 * std::log10(std::sqrt(err_with_dpd / ref_power));
    }

    return res;
}

QPair<double, double> MetricsEval::computeACPR(const std::vector<double>& freq,
                                               const std::vector<double>& psd,
                                               double BB,
                                               double deltaf,
                                               const QString& SigType)
{
    Q_UNUSED(SigType);

    ACPRResult full = computeACPR_full(freq, psd, BB, deltaf);
    return {full.lower_dB, full.upper_dB};
}

ACPRResult MetricsEval::computeACPR_full(const std::vector<double>& freq,
                                         const std::vector<double>& psd,
                                         double BB,
                                         double deltaf)
{
    ACPRResult res{};

    res.lower_dB = 0.0;
    res.upper_dB = 0.0;
    res.P_main_dBm = -std::numeric_limits<double>::infinity();
    res.P_lower_dBm = -std::numeric_limits<double>::infinity();
    res.P_upper_dBm = -std::numeric_limits<double>::infinity();
    res.P_main_mW = 0.0;
    res.P_lower_mW = 0.0;
    res.P_upper_mW = 0.0;

    if (freq.size() < 2 || freq.size() != psd.size())
        return res;

    // freq хранится в MHz.
    // BB и deltaf приходят в Hz.
    const double BB_MHz = BB / 1e6;
    const double df_MHz = deltaf / 1e6;

    const QPair<double, double> main_ch = {
        -BB_MHz / 2.0,
        BB_MHz / 2.0
    };

    const QPair<double, double> upper_ch = {
        -BB_MHz / 2.0 + df_MHz,
        BB_MHz / 2.0 + df_MHz
    };

    const QPair<double, double> lower_ch = {
        -BB_MHz / 2.0 - df_MHz,
        BB_MHz / 2.0 - df_MHz
    };

    for (int i = 0; i < static_cast<int>(freq.size()); ++i)
    {
        const double f = freq[i];

        // ВАЖНО:
        // psd[i] сейчас dBm/bin, то есть мощность конкретного FFT/Welch-бина.
        // Поэтому переводим dBm -> mW и просто суммируем бины.
        const double P_bin_mW = std::pow(10.0, psd[i] / 10.0);

        if (f >= lower_ch.first && f < lower_ch.second)
            res.P_lower_mW += P_bin_mW;

        if (f >= main_ch.first && f < main_ch.second)
            res.P_main_mW += P_bin_mW;

        if (f >= upper_ch.first && f < upper_ch.second)
            res.P_upper_mW += P_bin_mW;
    }

    auto mW_to_dBm_local = [](double p_mW) -> double
    {
        if (p_mW <= 0.0)
            return -std::numeric_limits<double>::infinity();

        return 10.0 * std::log10(p_mW);
    };

    res.P_main_dBm  = mW_to_dBm_local(res.P_main_mW);
    res.P_lower_dBm = mW_to_dBm_local(res.P_lower_mW);
    res.P_upper_dBm = mW_to_dBm_local(res.P_upper_mW);

    res.lower_dB =
        (res.P_main_mW > 0.0 && res.P_lower_mW > 0.0)
            ? 10.0 * std::log10(res.P_lower_mW / res.P_main_mW)
            : -std::numeric_limits<double>::infinity();

    res.upper_dB =
        (res.P_main_mW > 0.0 && res.P_upper_mW > 0.0)
            ? 10.0 * std::log10(res.P_upper_mW / res.P_main_mW)
            : -std::numeric_limits<double>::infinity();

    return res;
}

double MetricsEval::compute_av_P(const std::vector<std::complex<double>>& tx)
{
    // После привязки PA к Pin_sat/Pout_sat считаем:
    // |x|^2 = мощность в mW.
    return avgPower_mW(tx);
}

double MetricsEval::compute_av_P_dBm(const std::vector<std::complex<double>>& tx)
{
    return mW_to_dBm(avgPower_mW(tx));
}

double MetricsEval::compute_av_P_G(double Pin_mW, double Pout_mW)
{
    return 10.0 * std::log10(Pout_mW / Pin_mW);
}

double MetricsEval::computeIBO_dB(const std::vector<std::complex<double>>& paInput,
                                  double Pin_sat_dBm)
{
    const double Pin_avg_dBm = compute_av_P_dBm(paInput);
    return Pin_sat_dBm - Pin_avg_dBm;
}

double MetricsEval::computeOBO_dB(const std::vector<std::complex<double>>& paOutput,
                                  double Pout_sat_dBm)
{
    const double Pout_avg_dBm = compute_av_P_dBm(paOutput);
    return Pout_sat_dBm - Pout_avg_dBm;
}

std::vector<double> MetricsEval::hamming(int N)
{
    std::vector<double> w(N);

    if (N <= 1)
    {
        if (N == 1)
            w[0] = 1.0;
        return w;
    }

    for (int n = 0; n < N; ++n)
        w[n] = 0.54 - 0.46 * std::cos(2.0 * M_PI * n / (N - 1));

    return w;
}

void MetricsEval::computePSDWelch(
    const std::vector<std::complex<double>>& tx,
    double Fs,
    double oversample,
    std::vector<double>& freq,
    std::vector<double>& psd)
{
    const int winSize = PSD_WIN_SIZE;
    const int overlap = PSD_OVERLAP;
    const int step = winSize - overlap;
    const int nfft = PSD_NFFT;

    if (tx.size() < static_cast<size_t>(winSize))
        return;

    std::vector<double> localWindow = hamming(winSize);

    const int segments = 1 + (static_cast<int>(tx.size()) - winSize) / step;

    psd = std::vector<double>(nfft, 0.0);

    std::vector<std::complex<double>> buffer(nfft);

    for (int k = 0; k < segments; ++k)
    {
        const int start = k * step;

        std::fill(buffer.begin(), buffer.end(), std::complex<double>(0.0, 0.0));

        for (int i = 0; i < winSize; ++i)
            buffer[i] = tx[start + i] * localWindow[i];

        fftw8192.fftInPlace(buffer);

        for (int i = 0; i < nfft; ++i)
            psd[i] += std::norm(buffer[i]);
    }

    for (int i = 0; i < nfft; ++i)
        psd[i] /= static_cast<double>(segments);

    double U = 0.0;

    for (double w : localWindow)
        U += w * w;

    U /= static_cast<double>(winSize);

    const double Fs_eff = Fs * oversample;

    // Если |x|^2 = mW, то после этой нормировки PSD = mW/Hz.
    const double scale = Fs_eff * winSize * U;

    for (int i = 0; i < nfft; ++i)
        psd[i] /= scale;

    std::rotate(psd.begin(), psd.begin() + nfft / 2, psd.end());

    freq.resize(nfft);

    for (int i = 0; i < nfft; ++i)
        freq[i] = (i - nfft / 2) * Fs_eff / nfft;

    const double eps = 1e-30;
    double df_Hz = Fs_eff / nfft;

    for (int i = 0; i < nfft; ++i)
    {
        freq[i] = freq[i] / 1e6;

        double p_bin_mW = psd[i] * df_Hz; // PSD[mW/Hz] * Hz = mW

        psd[i] = 10.0 * std::log10(p_bin_mW + eps); // dBm/bin
    }
}

double MetricsEval::computePAPR_dB(const std::vector<std::complex<double>>& x)
{
    if (x.empty())
        return 0.0;

    double pavg = 0.0;
    double ppeak = 0.0;

    for (const auto& s : x)
    {
        double p = std::norm(s);
        pavg += p;

        if (p > ppeak)
            ppeak = p;
    }

    pavg /= static_cast<double>(x.size());

    if (pavg <= 0.0)
        return 0.0;

    return 10.0 * std::log10(ppeak / pavg);
}
