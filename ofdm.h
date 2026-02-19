#ifndef OFDM_H
#define OFDM_H
#include <QVector>
#include <complex>
#include <cmath>
#include <QRandomGenerator>
#include <fft.h>

struct OfdmParams {
    int Nfft = 32;
    double fs = 0.0;
    double fc = 0.0;
    int GB_DC = 0;
    int GB_Nyq = 0;
    int CP = 0;
    double SNR_dB = 100.0;
    int oversampling = 1;
};

struct OfdmResult {
    QVector<std::complex<double>> tx;                     // выходной сигнал
    QVector<double> t;                      // время
    QVector<std::complex<double>> currentNoise;
    double fc;
    double BB = 0.0;                        // полоса
    double Tsym = 0.0;                      // длительность символа
    int Nactive = 0;                        // число активных поднесущих
};

class OFDM
{
public:
    OFDM();
    OfdmResult makeOfdm(const QVector<std::complex<double>>& symbols,
                        const OfdmParams& p);
    QVector<std::complex<double>> ofdm_demodulate(
        const QVector<std::complex<double>> &rx,
        const QVector<std::complex<double>> &sym_clean,
        const OfdmParams& p);
    void changeFc(OfdmResult &x, OfdmParams& p);
    void changeAwgn(OfdmResult &x, OfdmParams& p);
private:
    QVector<QVector<std::complex<double>>> ofdmSubcarrierMapping(
        const QVector<std::complex<double>>& dataSymbols,
        const int Nfft, const int GB_DC, const int GB_Nyq, int& Nactive);
    QVector<std::complex<double>> ofdm_subcarrier_demapping(
        const QVector<QVector<std::complex<double>>>& X,
        const int Nfft, const int GB_DC, const int GB_Nyq, const int NumSym);
    void addAwgn(OfdmResult &x, double SNR_dB);
};

#endif // OFDM_H
