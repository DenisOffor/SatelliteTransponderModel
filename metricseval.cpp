#include "metricseval.h"

MetricsEval::MetricsEval(): fftw8192(FFT_SIZE) {
    window = hamming(WINDOW_SISE);
}

std::vector<std::complex<double>> MetricsEval::normalizeSignal(
    const std::vector<std::complex<double>>& tx)
{
    if (tx.empty())
        return {};

    double sum = 0.0;
    for (const auto& s : tx)
        sum += std::norm(s);   // |x|^2

    double rms = std::sqrt(sum / tx.size());

    if (rms == 0.0)
        return tx;

    std::vector<std::complex<double>> tx_norm(tx.size());
    for (int i = 0; i < tx.size(); ++i)
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

    for (size_t i = 0; i < ref.size(); ++i) {
        num += std::conj(ref[i]) * test[i];
        den += std::norm(ref[i]);
    }

    if (den <= 0.0)
        return 0.0;

    std::complex<double> a = num / den;

    double errPow = 0.0;
    double refPow = 0.0;

    for (size_t i = 0; i < ref.size(); ++i) {
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
    std::vector<double>& psd_tx,
    std::vector<std::complex<double>>& sig_buff)
{
    if(tx.size() < WINDOW_SISE) {
        if(sig_buff.size() > 4 * WINDOW_SISE)
            sig_buff.clear();

        for(int i = 0; i < tx.size(); ++i)
            sig_buff.push_back(tx[i]);

        if(sig_buff.size() < 4 * WINDOW_SISE)
            return;
        else {
            freq.clear();
            psd_tx.clear();
            auto tx_norm = normalizeSignal(sig_buff);
            computePSDWelch(tx_norm, Fs, oversample, freq, psd_tx);
            return;
        }
    }


    freq.clear();
    psd_tx.clear();
    auto tx_norm = normalizeSignal(tx);
    computePSDWelch(tx_norm, Fs, oversample, freq, psd_tx);
}

QPair<double, double> MetricsEval::Calc_BER(std::vector<Symbols> &symbols)
{
    double errors_noDPD = 0;
    double errors_withDPD = 0;
    for(int n = 0; n < symbols.size(); ++n) {
        for(int i = 0; i < symbols[n].data_rx.size(); ++i) {
            if(symbols[n].data_rx[i] != symbols[n].data_tx[i])
                errors_noDPD++;
            if(symbols[n].data_rx_with_DPD[i] != symbols[n].data_tx[i])
                errors_withDPD++;
        }
    }
    return QPair<double, double> (1.0 * errors_noDPD / (symbols.size() * symbols[0].data_tx.size()), 1.0 * errors_withDPD / (symbols.size() * symbols[0].data_tx.size()));
}

QPair<double, double> MetricsEval::Calc_EVM(const std::vector<Symbols>& frames)
{
    double err_no_dpd = 0.0;
    double err_with_dpd = 0.0;
    double ref_power = 0.0;

    size_t N = 0;

    for(const auto& f : frames)
    {
        size_t n = std::min({f.tr_sym_clean.size(),
                             f.rec_sym_noisy.size(),
                             f.rec_sym_noisy_with_DPD.size()});

        for(size_t i = 0; i < n; ++i)
        {
            const auto& ref = f.tr_sym_clean[i];

            auto err1 = f.rec_sym_noisy[i] - ref;
            auto err2 = f.rec_sym_noisy_with_DPD[i] - ref;

            err_no_dpd += std::norm(err1);
            err_with_dpd += std::norm(err2);

            ref_power += std::norm(ref);

            ++N;
        }
    }

    QPair<double, double> res;

    if(ref_power > 0)
    {
        res.first = 20*std::log10(std::sqrt(err_no_dpd / ref_power));
        res.second = 20*std::log10(std::sqrt(err_with_dpd / ref_power));
    }

    return res;
}

QPair<double, double> MetricsEval::computeACPR(const std::vector<double> &freq, const std::vector<double> &psd,
                                               double BB, double deltaf, const QString &SigType)
{
    if (SigType == "FDMA")
        return {0.0, 0.0};

    if (freq.size() < 2 || freq.size() != psd.size())
        return {0.0, 0.0};

    const double f_step = std::abs(freq[1] - freq[0]);

    const QPair<double, double> main_ch  = {-BB / (1e6 * 2.0),  BB / (1e6 * 2.0)};
    const QPair<double, double> upper_ch = {-BB / (1e6 * 2.0) + deltaf * 1.3 / 1e6,  BB / (1e6 * 2.0) + deltaf * 1.3 / 1e6};
    const QPair<double, double> lower_ch = {-BB / (1e6 * 2.0) - deltaf * 1.3 / 1e6,  BB / (1e6 * 2.0) - deltaf * 1.3 / 1e6};

    double P_lower_ch = 0.0;
    double P_upper_ch = 0.0;
    double P_main_ch  = 0.0;

    for (int i = 0; i < static_cast<int>(freq.size()); ++i) {
        const double f = freq[i];
        double psd_lin = std::pow(10.0, psd[i] / 10.0);

        if (f >= lower_ch.first && f <= lower_ch.second)
            P_lower_ch += psd_lin * f_step;
        else if (f >= main_ch.first && f <= main_ch.second)
            P_main_ch += psd_lin * f_step;
        else if (f >= upper_ch.first && f <= upper_ch.second)
            P_upper_ch += psd_lin * f_step;
    }

    if (P_main_ch <= 0.0)
        return {0.0, 0.0};

    const double ACPR_lower = (P_lower_ch > 0.0)
                                  ? 10.0 * std::log10(P_lower_ch / P_main_ch)
                                  : -std::numeric_limits<double>::infinity();

    const double ACPR_upper = (P_upper_ch > 0.0)
                                  ? 10.0 * std::log10(P_upper_ch / P_main_ch)
                                  : -std::numeric_limits<double>::infinity();

    return {ACPR_lower, ACPR_upper};
}

double MetricsEval::compute_av_P(const std::vector<std::complex<double> > &tx)
{
    if (tx.empty()) return 0.0;

    double sum = 0.0;
    for (const auto& s : tx)
        sum += std::norm(s);   // |x|^2

    return sum / (tx.size());
}

double MetricsEval::compute_av_P_G(double Pin, double Pout)
{
    if (Pin <= 0.0 || Pout <= 0.0)
        return -std::numeric_limits<double>::infinity();

    return 10.0 * std::log10(Pout / Pin);
}

std::vector<double> MetricsEval::hamming(int N)
{
    std::vector<double> w(N);
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
    const int winSize = 2048;
    const int overlap = 1024;
    const int step = winSize - overlap;
    const int nfft = 8192;

    if (tx.size() < winSize)
        return;


    // MATLAB-совместимое число сегментов
    int segments = 1 + (tx.size() - winSize) / step;

    psd = std::vector<double>(nfft, 0.0);

    std::vector<std::complex<double>> buffer(nfft);
    //QVector<std::complex<double>> spectrum(nfft);

    for (int k = 0; k < segments; ++k)
    {
        int start = k * step;

        std::fill(buffer.begin(), buffer.end(), 0.0);

        for (int i = 0; i < winSize; ++i)
            buffer[i] = tx[start + i] * window[i];

        fftw8192.fftInPlace(buffer);

        for (int i = 0; i < nfft; ++i)
            psd[i] += std::norm(buffer[i]);
    }

    // усреднение по сегментам
    for (int i = 0; i < nfft; ++i)
        psd[i] /= segments;

    // мощность окна
    double U = 0.0;
    for (double w : window)
        U += w * w;

    U /= winSize;   // MATLAB именно так делает

    // НОРМИРОВКА (ключевой момент!)
    double scale = Fs * winSize * U;

    for (int i = 0; i < nfft; ++i)
        psd[i] /= scale;

    // fftshift
    std::rotate(psd.begin(), psd.begin() + nfft/2, psd.end());

    // Частотная ось
    freq.resize(nfft);
    for (int i = 0; i < nfft; ++i)
        freq[i] = (i - nfft/2) * Fs * oversample / nfft;

    // перевод в dB
    const double eps = 1e-20;
    for (int i = 0; i < nfft; ++i) {
        freq[i] = freq[i] / 1e6;
        psd[i] = 10.0 * std::log10(psd[i] + eps);
    }
}

double MetricsEval::computePAPR_dB(const std::vector<std::complex<double>>& x)
{
    if (x.empty()) return 0.0;

    double pavg = 0.0;
    double ppeak = 0.0;

    for (const auto& s : x) {
        double p = std::norm(s);
        pavg += p;
        if (p > ppeak) ppeak = p;
    }

    pavg /= x.size();
    if (pavg <= 0.0) return 0.0;

    return 10.0 * std::log10(ppeak / pavg);
}
