#ifndef FDMA_H
#define FDMA_H
#include <complex>
#include <cmath>
#include <QRandomGenerator>
#include <fft.h>
#include "sc.h"
#include "HelpfullStructs.h"

struct FdmaParams
{
    int FDMA_f_carrier;        // стартовая
    int FDMA_symrate;
    int FDMA_num_subcarriers;
    int FDMA_step_carrier;
    int oversampling;
    int fs;
    int SNRSig;

    double rolloff = 0.35;
    int filterSpan = 6;
    QString filterType = "rrc";
};

struct FdmaResult
{
    std::vector<std::complex<double>> tx;
    std::vector<double> t;
    std::vector<std::complex<double>> currentNoise;
    double totalBandwidth;

    std::vector<ScResult> perCarrierResults;
};

class FDMA
{
public:
    FDMA(SC& scGenerator);

    FdmaResult generate(
        const std::vector<Symbols> symbolsPerCarrier,
        const FdmaParams& p, const ScParams& sc_p);
    std::vector<std::vector<std::complex<double>>> demodulate(
    const std::vector<std::complex<double>>& rxSignal, std::vector<Symbols>& symbols,
        const FdmaParams& p, const ScParams& sc_p);

private:
    SC& sc;
};

#endif // FDMA_H
