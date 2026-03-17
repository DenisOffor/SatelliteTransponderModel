#ifndef DPD_H
#define DPD_H

#include <vector>
#include <complex>
#include <Eigen/Dense>

class DPD {
public:
    void train(const std::vector<std::complex<double>> &pa_input,
               const std::vector<std::complex<double>> &pa_output,
               int P, int M);

    std::vector<std::complex<double>> applyPreDistortion(
        const std::vector<std::complex<double>> &input_signal,
        int P, int M);

private:
    std::vector<std::complex<double>> coeffs; // коэффициенты DPD
    double scale = 1.0;                       // внутренний масштаб
    std::vector<std::complex<double>> normalizeSignal(const std::vector<std::complex<double>> &sig);
};

#endif // DPD_H
