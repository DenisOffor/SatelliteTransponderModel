#ifndef DPD_H
#define DPD_H

#include <vector>
#include <complex>
#include <Eigen/Dense>
#include <QDebug>
#include "HelpfullStructs.h"
#include "QElapsedTimer"

using namespace Eigen;

class DPD {
public:
    std::vector<std::complex<double>> applyMP(const std::vector<std::complex<double>>& sig, const Source& source);
    std::vector<std::complex<double>> applyGMP(const std::vector<std::complex<double>>& sig, const Source& source);
    std::vector<std::complex<double>> coeffs;
    void train(const std::vector<std::complex<double>>& pa_input, const std::vector<std::complex<double>>& pa_output,
               const Source& source);
private:
    double Gpeak;
    double Grms;
    double computePeakGain(const std::vector<std::complex<double>>& pa_input,
                                const std::vector<std::complex<double>>& pa_output);
    std::vector<std::complex<double>> normalizeByGain(const std::vector<std::complex<double>>& pa_output,
        double Gpeak);
    double computeGrms(const std::vector<std::complex<double>>& pa_input,
                            const std::vector<std::complex<double>>& pa_output);
    double computeRMS(const std::vector<std::complex<double>>& sig);
    VectorXcd DPDsolve_least_squares(const MatrixXcd& A, const VectorXcd& b);
    MatrixXcd make_MP_mat(const std::vector<std::complex<double>>& x, const int P, const int M);
    MatrixXcd make_GMP_mat(const std::vector<std::complex<double>>& x, const int P, const int M, const int L_lag, const int L_lead);
    VectorXcd make_goal(const  std::vector<std::complex<double>>& y, const Source& source);
};

#endif // DPD_H
