#ifndef FDMA_H
#define FDMA_H
#include <QVector>
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
    QVector<std::complex<double>> tx;
    QVector<double> t;
    double totalBandwidth;

    QVector<ScResult> perCarrierResults;
};

class FDMA
{
public:
    FDMA(SC& scGenerator);

    FdmaResult generate(
        const QVector<Symbols> symbolsPerCarrier,
        const FdmaParams& p);

private:
    SC& sc;
};

#endif // FDMA_H
