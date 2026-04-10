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

// void MetricsEval::comparePSD(
//     const std::vector<std::complex<double>>& tx,
//     const std::vector<std::complex<double>>& rx,
//     double Fs,
//     double oversample,
//     std::vector<double>& freq,
//     std::vector<double>& psd_tx,
//     std::vector<double>& psd_rx)
// {
//     if(sig_buff.size() > WINDOW_SISE * 2) {
//         sig_buff.clear();
//     }

//     for(int i = 0; i < tx.size())


//     freq.clear();
//     psd_rx.clear();
//     psd_tx.clear();
//     auto tx_norm = normalizeSignal(tx);
//     auto rx_norm = normalizeSignal(rx);

//     computePSDWelch(tx_norm, Fs, oversample, freq, psd_tx);
//     computePSDWelch(rx_norm, Fs, oversample, freq, psd_rx);
// }

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
