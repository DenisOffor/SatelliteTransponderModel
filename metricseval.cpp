#include "metricseval.h"

MetricsEval::MetricsEval() {}

QVector<std::complex<double>> MetricsEval::normalizeSignal(
    const QVector<std::complex<double>>& tx)
{
    if (tx.isEmpty())
        return {};

    double sum = 0.0;
    for (const auto& s : tx)
        sum += std::norm(s);   // |x|^2

    double rms = std::sqrt(sum / tx.size());

    if (rms == 0.0)
        return tx;

    QVector<std::complex<double>> tx_norm(tx.size());
    for (int i = 0; i < tx.size(); ++i)
        tx_norm[i] = tx[i] / rms;

    return tx_norm;
}

void MetricsEval::comparePSD(
    const QVector<std::complex<double>>& tx,
    const QVector<std::complex<double>>& rx,
    double Fs,
    QVector<double>& freq,
    QVector<double>& psd_tx,
    QVector<double>& psd_rx)
{
    auto tx_norm = normalizeSignal(tx);
    auto rx_norm = normalizeSignal(rx);

    computePSDWelch(tx_norm, Fs, freq, psd_tx);
    computePSDWelch(rx_norm, Fs, freq, psd_rx);
}

QVector<double> MetricsEval::hamming(int N)
{
    QVector<double> w(N);
    for (int n = 0; n < N; ++n)
        w[n] = 0.54 - 0.46 * std::cos(2.0 * M_PI * n / (N - 1));
    return w;
}

void MetricsEval::computePSDWelch(
    const QVector<std::complex<double>>& tx,
    double Fs,
    QVector<double>& freq,
    QVector<double>& psd)
{
    const int winSize = 2048;
    const int overlap = 1024;
    const int step = winSize - overlap;
    const int nfft = 8192;

    if (tx.size() < winSize)
        return;

    auto window = hamming(winSize);

    // MATLAB-совместимое число сегментов
    int segments = 1 + (tx.size() - winSize) / step;

    psd = QVector<double>(nfft, 0.0);

    FFT myfft(nfft);
    QVector<std::complex<double>> buffer(nfft);

    for (int k = 0; k < segments; ++k)
    {
        int start = k * step;

        std::fill(buffer.begin(), buffer.end(), 0.0);

        for (int i = 0; i < winSize; ++i)
            buffer[i] = tx[start + i] * window[i];

        auto spectrum = myfft.fft(buffer);

        for (int i = 0; i < nfft; ++i)
            psd[i] += std::norm(spectrum[i]);
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
        freq[i] = (i - nfft/2) * Fs / nfft;

    // перевод в dB
    const double eps = 1e-20;
    for (int i = 0; i < nfft; ++i) {
        freq[i] = freq[i] / 1e6;
        psd[i] = 10.0 * std::log10(psd[i] + eps);
    }
}
