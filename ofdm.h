#ifndef OFDM_H
#define OFDM_H
#include <complex>
#include <cmath>
#include <QRandomGenerator>
#include <fft.h>

struct OfdmParams {
    int Nfft = 32;
    double fs = 0.0;
    int GB_DC = 0;
    int GB_Nyq = 0;
    int CP = 0;
    double SNR_dB = 100.0;
    int oversampling = 1;
};

struct OfdmResult {
    std::vector<std::complex<double>> tx;                     // выходной сигнал
    std::vector<double> t;                      // время

    double fc;
    double BB = 0.0;                        // полоса
    double Tsym = 0.0;                      // длительность символа
    int Nactive = 0;                        // число активных поднесущих
};

class OFDM
{
public:
    OFDM();
    OfdmResult makeOfdm(const std::vector<std::complex<double>>& symbols,
                        const OfdmParams& p);
    std::vector<std::complex<double>> ofdm_demodulate(
        const std::vector<std::complex<double>> &rx,
        const std::vector<std::complex<double>> &sym_clean,
        const OfdmParams& p);
private:
    std::vector<std::vector<std::complex<double>>> ofdmSubcarrierMapping(
        const std::vector<std::complex<double>>& dataSymbols,
        const int Nfft, const int GB_DC, const int GB_Nyq, int& Nactive);
    std::vector<std::complex<double>> ofdm_subcarrier_demapping(
        const std::vector<std::vector<std::complex<double>>>& X,
        const int Nfft, const int GB_DC, const int GB_Nyq, const int NumSym);
};

#endif // OFDM_H
