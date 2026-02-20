#ifndef METRICSEVAL_H
#define METRICSEVAL_H
#include "HelpfullStructs.h"

class MetricsEval
{
public:
    MetricsEval();
    void computePSDWelch(const QVector<std::complex<double>>& tx, double Fs,
                         QVector<double>& freq, QVector<double>& psd);
    QVector<std::complex<double>> normalizeSignal(const QVector<std::complex<double>>& tx);
    QVector<double> hamming(int N);
    void comparePSD(const QVector<std::complex<double>>& tx,
                    const QVector<std::complex<double>>& rx, double Fs,
                    QVector<double>& freq, QVector<double>& psd_tx, QVector<double>& psd_rx);
};

#endif // METRICSEVAL_H
