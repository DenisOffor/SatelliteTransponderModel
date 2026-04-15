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
    void computePSD(
        const std::vector<std::complex<double>>& tx,
        double Fs,
        double oversample,
        std::vector<double>& freq,
        std::vector<double>& psd_tx,
        std::vector<std::complex<double>>& sig_buff);
    QPair<double, double> Calc_BER(std::vector<Symbols>& symbols);
    QPair<double, double> Calc_EVM(const std::vector<Symbols>& frames);
    QPair<double, double> computeACPR(const std::vector<double> &freq, const std::vector<double> &psd, double BB,
                                                   double deltaf, const QString &SigType);
    double compute_av_P(const std::vector<std::complex<double>>& tx);
    double compute_av_P_G(double Pin, double Pout);
    double computePAPR_dB(const std::vector<std::complex<double>>& x);
    double computeNMSE_dB(const std::vector<std::complex<double>>& ref,
                                const std::vector<std::complex<double>>& test);
};


#endif // METRICSEVAL_H
