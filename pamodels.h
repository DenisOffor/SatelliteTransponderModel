#ifndef PAMODELS_H
#define PAMODELS_H
#include "HelpfullStructs.h"

class PAModels
{
public:
    PAModels();
    void SalehModel(QVector<std::complex<double>>& sig, QVector<double>& Coeffs, int& linear_gain_dB, int& IBO_dB);
    void RappModel(QVector<std::complex<double>>& sig, QVector<double>& Coeffs, int& linear_gain_dB, int& IBO_dB);
    void GhorbaniModel(QVector<std::complex<double>>& sig, QVector<double>& Coeffs, int& linear_gain_dB, int& IBO_dB);
    void MemoryPolinomalModel(QVector<std::complex<double>>& sig, QVector<double>& Coeffs, int& linear_gain_dB, int& IBO_dB);
    void apply_IBO(QVector<std::complex<double>>& tx, double IBO_dB);
};

#endif // PAMODELS_H
