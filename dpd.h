#ifndef DPD_H
#define DPD_H

#include <complex>
#include <cmath>
#include <Eigen/Dense>

class DPD
{
public:

    DPD();

    void train(const std::vector<std::complex<double>>& pa_input,
               const std::vector<std::complex<double>>& pa_output);

    std::vector<std::complex<double>> predistort(
        const std::vector<std::complex<double>>& input);

    void setOrder(int order);
    void setMemory(int memory);

private:

    int K; // nonlinear order
    int M; // memory depth

    std::vector<std::complex<double>> coeffs;

    std::complex<double> basis(
        const std::complex<double>& x,
        int k);
};

#endif // DPD_H
