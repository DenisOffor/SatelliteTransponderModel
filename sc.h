#ifndef SC_H
#define SC_H
#include <QVector>
#include <complex>
#include <cmath>
#include <QRandomGenerator>
#include <fft.h>

struct ScParams {
    int SC_f_carrier = 0;        // можно не использовать отдельно
    int SC_symrate = 0;          // Rs
    double SC_rolloff = 0.35;
    int SC_filter_length = 6;    // span
    QString SC_FilterType = "RRC";

    double fs = 0.0;
    double fc = 0.0;
    double SNR_dB = 100.0;
    int oversampling = 1;
};

struct ScResult {
    QVector<std::complex<double>> tx;
    QVector<double> t;
    QVector<std::complex<double>> currentNoise;

    double bandwidth = 0.0;
    double Tsym = 0.0;
    double fc = 0.0;
};

class SC
{
public:
    SC();

    ScResult makeSc(const QVector<std::complex<double>>& symbols,
                    const ScParams& p);
    QVector<std::complex<double>> demodulateSignal(
        const QVector<std::complex<double>>& tx_signal,
        const ScParams& p,
        ScParams& updatedParams);
    void changeAwgn(ScResult &x, ScParams &p);

private:
    QVector<std::complex<double>> upsample(
        const QVector<std::complex<double>>& in, int sps);

    QVector<double> rrcFilter(int span, int sps, double beta);

    QVector<std::complex<double>> filterSignal(
        const QVector<std::complex<double>>& x,
        const QVector<double>& h);
    QVector<std::complex<double>> polyphaseFilter(
        const QVector<std::complex<double>>& symbols,
        const QVector<double>& h,
        int sps);
    QVector<std::complex<double>> matchedFilterDecimate(
        const QVector<std::complex<double>>& x,
        const QVector<double>& h,
        int sps);

    void addAwgn(ScResult &x, double SNR_dB);
};

#endif // SC_H
