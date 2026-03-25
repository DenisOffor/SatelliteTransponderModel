#ifndef DPD_H
#define DPD_H

#include <vector>
#include <complex>
#include <Eigen/Dense>
#include <QDebug>

using namespace Eigen;

class DPD {
public:
    std::vector<std::complex<double>> applyMP(const std::vector<std::complex<double>>& sig, const int P, const int M);
    std::vector<std::complex<double>> coeffs;
    void train(const std::vector<std::complex<double>>& pa_input, const std::vector<std::complex<double>>& pa_output,
               const int P, const int M);
private:
    double computePeakGain(const std::vector<std::complex<double>>& pa_input,
                                const std::vector<std::complex<double>>& pa_output);
    std::vector<std::complex<double>> normalizeByPeakGain(
        const std::vector<std::complex<double>>& pa_output,
        double Gpeak);
    VectorXcd DPDsolve_least_squares(const MatrixXcd& A, const VectorXcd& b);
    MatrixXcd make_MP_mat(const std::vector<std::complex<double>>& x, const int P, const int M);
    VectorXcd make_goal(const  std::vector<std::complex<double>>& y, const int M);
};

#endif // DPD_H
