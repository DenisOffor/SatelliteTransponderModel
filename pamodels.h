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
    void HammersteinModel(std::vector<std::complex<double>>& sig, QString Static_model, std::vector<double>& Coeffs,
                     std::vector<double>& FIR_Coeffs, int& linear_gain_dB, int& IBO_dB);
    void WHModel(std::vector<std::complex<double>>& sig, QString Static_model, std::vector<double>& Coeffs,
                     std::vector<double>& FIR_Coeffs, int& linear_gain_dB, int& IBO_dB);
    double find_Asat_Ghorbani(const std::vector<double>& c, double gain_linear);
    void ApplyFIRWithMemory(std::vector<std::complex<double>>& signal, double alpha, int numTaps);
    void ScaleToRMS_forPA(Source& source, GlobalResults& CurRes);
    void scaleToRMS(std::vector<std::complex<double>>& x, double target_rms);
};

#endif // PAMODELS_H
