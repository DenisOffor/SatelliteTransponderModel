#ifndef PAMODELS_H
#define PAMODELS_H
#include "HelpfullStructs.h"

class PAModels
{
public:
    PAModels();
    void SalehModel(std::vector<std::complex<double>>& sig, std::vector<double>& Coeffs, int& linear_gain_dB, int& IBO_dB);
    void RappModel(std::vector<std::complex<double>>& sig, std::vector<double>& Coeffs, int& linear_gain_dB, int& IBO_dB);
    void GhorbaniModel(std::vector<std::complex<double>>& sig, std::vector<double>& Coeffs, int& linear_gain_dB, int& IBO_dB);
    void WienerModel(std::vector<std::complex<double>>& sig, QString Static_model, std::vector<double>& Coeffs,
                     std::vector<double>& FIR_Coeffs, int& linear_gain_dB, int& IBO_dB);
    void apply_IBO(std::vector<std::complex<double>>& tx, double IBO_dB, double A_sat);
    double find_Asat_Ghorbani(const std::vector<double>& c, double gain_linear);
    void ApplyFIRWithMemory(std::vector<double>& amplitude,std::vector<double>& phase,
                            const std::vector<double>& FIR_Coefs, int numTaps);
};

#endif // PAMODELS_H
