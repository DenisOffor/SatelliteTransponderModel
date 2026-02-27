#ifndef METRICSEVAL_H
#define METRICSEVAL_H
#include "HelpfullStructs.h"

#define FFT_SIZE 8192
#define WINDOW_SISE 2048

class MetricsEval
{
private:
    FFT fftw8192;
    std::vector<double> window;
public:
    MetricsEval();
    void computePSDWelch(const std::vector<std::complex<double>>& tx, double Fs, double oversample,
                         std::vector<double>& freq, std::vector<double>& psd);
    std::vector<std::complex<double>> normalizeSignal(const std::vector<std::complex<double>>& tx);
    std::vector<double> hamming(int N);
    void comparePSD(const std::vector<std::complex<double>>& tx,
                    const std::vector<std::complex<double>>& rx, double Fs, double oversample,
                    std::vector<double>& freq, std::vector<double>& psd_tx, std::vector<double>& psd_rx);
};


#endif // METRICSEVAL_H
