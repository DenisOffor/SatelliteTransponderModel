#ifndef METRICSEVAL_H
#define METRICSEVAL_H
#include "HelpfullStructs.h"

#define FFT_SIZE 8192
#define WINDOW_SISE 2048

struct ACPRResult
{
    double lower_dB = 0.0;
    double upper_dB = 0.0;

    double P_main_dBm  = 0.0;
    double P_lower_dBm = 0.0;
    double P_upper_dBm = 0.0;

    double P_main_mW  = 0.0;
    double P_lower_mW = 0.0;
    double P_upper_mW = 0.0;
};

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

    double compute_av_P_dBm(const std::vector<std::complex<double>>& tx);
    double computeIBO_dB(const std::vector<std::complex<double>>& paInput, double Pin_sat_dBm);
    double computeOBO_dB(const std::vector<std::complex<double>>& paOutput, double Pout_sat_dBm);
    ACPRResult computeACPR_full(const std::vector<double>& freq, const std::vector<double>& psd, double BB, double deltaf);
};


#endif // METRICSEVAL_H
