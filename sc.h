#ifndef SC_H
#define SC_H

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
    std::vector<std::complex<double>> tx;
    std::vector<double> t;

    double bandwidth = 0.0;
    double Tsym = 0.0;
    double fc = 0.0;
};

class SC
{
public:
    SC();

    ScResult makeSc(const std::vector<std::complex<double>>& symbols,
                    const ScParams& p);
    std::vector<std::complex<double>> demodulateSignal(
        const std::vector<std::complex<double>>& tx_signal, const std::vector<std::complex<double>> &sym_clean,
        const ScParams& p,
        ScParams& updatedParams);

private:
    std::vector<std::complex<double>> upsample(
        const std::vector<std::complex<double>>& in, int sps);

    std::vector<double> rrcFilter(int span, int sps, double beta);

    std::vector<std::complex<double>> filterSignal(
        const std::vector<std::complex<double>>& x,
        const std::vector<double>& h);
    std::vector<std::complex<double>> polyphaseFilter(
        const std::vector<std::complex<double>>& symbols,
        const std::vector<double>& h,
        int sps);
    std::vector<std::complex<double>> matchedFilterDecimate(
        const std::vector<std::complex<double>>& x,
        const std::vector<double>& h,
        int sps);
};

#endif // SC_H
