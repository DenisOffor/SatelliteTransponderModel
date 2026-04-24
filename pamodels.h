#ifndef PAMODELS_H
#define PAMODELS_H
#include "HelpfullStructs.h"

class PAModels
{
public:
    PAModels();

    void SalehModel(std::vector<std::complex<double>>& sig, Source& source);
    void RappModel(std::vector<std::complex<double>>& sig, Source& source);
    void GhorbaniModel(std::vector<std::complex<double>>& sig, Source& source);
    void WienerModel(std::vector<std::complex<double>>& sig, QString Static_model, std::vector<double>& FIR_Coeffs, Source& source);
    void HammersteinModel(std::vector<std::complex<double>>& sig, QString Static_model, std::vector<double>& FIR_Coeffs, Source& source);
    void WHModel(std::vector<std::complex<double>>& sig, QString Static_model, std::vector<double>& FIR_Coeffs, Source& source);
    void ApplyPA(std::vector<std::complex<double>>& sig, Source& source);

    double find_Asat_Ghorbani(const std::vector<double>& c, double gain_linear);
    void ApplyFIRWithMemory(std::vector<std::complex<double>>& signal, double alpha, int numTaps);
    void ScaleToRMS_forPA(std::vector<std::complex<double>>& sig, Source& source);
    void scaleToRMS(std::vector<std::complex<double>>& x, double target_rms);
};

#endif // PAMODELS_H
