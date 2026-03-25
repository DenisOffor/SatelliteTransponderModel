#include "dpd.h"
#include <cmath>

void DPD::train(const std::vector<std::complex<double>> &pa_input,
                const std::vector<std::complex<double> > &pa_output, const int P, const int M)
{
    double Gpeak = computePeakGain(pa_input, pa_output);
    std::vector<std::complex<double>> pa_output_norm =
        normalizeByPeakGain(pa_output, Gpeak);

    VectorXcd a = DPDsolve_least_squares(
        make_MP_mat(pa_output_norm, P, M),
        make_goal(pa_input, M)
        );

    coeffs.resize(M * (P + 1) / 2);
    for (int i = 0; i < M * ((P + 1) / 2); ++i)
        coeffs[i] = a(i);
}

double DPD::computePeakGain(const std::vector<std::complex<double>>& pa_input,
                            const std::vector<std::complex<double>>& pa_output)
{
    double max_in = 0.0;
    double max_out = 0.0;

    for (const auto& s : pa_input)
        max_in = std::max(max_in, std::abs(s));

    for (const auto& s : pa_output)
        max_out = std::max(max_out, std::abs(s));

    if (max_in == 0.0)
        return 1.0;

    return max_out / max_in;
}

std::vector<std::complex<double>> DPD::normalizeByPeakGain(
    const std::vector<std::complex<double>>& pa_output,
    double Gpeak)
{
    std::vector<std::complex<double>> y_norm(pa_output.size());

    if (Gpeak == 0.0)
        return y_norm;

    for (size_t i = 0; i < pa_output.size(); ++i)
        y_norm[i] = pa_output[i] / Gpeak;

    return y_norm;
}

MatrixXcd DPD::make_MP_mat(const std::vector<std::complex<double>>& x, const int P, const int M)
{
    // y = Phi(x...)*a
    // исключили первые M-1 отсчетов, чтобы не уходить в x[-1] и тд, то есть отрицательные отсчеты
    int N = x.size();
    int rows = N - M + 1;
    int cols = M * (P + 1) / 2;
    MatrixXcd Phi_LS(rows, cols);

    for(int n = M - 1; n < N; ++n) {
        size_t cur_row = n - M + 1;
        size_t col_idx = 0;

        for(int p = 1; p <= P; p += 2)
            for(int m = 0; m < M; ++m)
                Phi_LS(cur_row, col_idx++) = x[n - m] * std::pow(std::abs(x[n - m]), p - 1);
    }

    return Phi_LS;
}

VectorXcd DPD::make_goal(const std::vector<std::complex<double>>& y, const int M) {
    int N = y.size();
    VectorXcd y_LS(N - M + 1);

    size_t num_el = 0;
    for(int n = M - 1; n < N; ++n)
        y_LS(num_el++) = y[n];

    return y_LS;
}

VectorXcd DPD::DPDsolve_least_squares(const MatrixXcd& Phi, const VectorXcd& goal) {
    // Нахождение решения МНК через нормальное уравнение
    return (Phi.adjoint() * Phi).ldlt().solve(Phi.adjoint() * goal);
}

std::vector<std::complex<double>> DPD::applyMP(
    const std::vector<std::complex<double>>& sig, int P, int M)
{
    std::vector<std::complex<double>> x_pre(sig.size(), std::complex<double>(0.0, 0.0));

    for (int i = 0; i < M - 1 && i < static_cast<int>(sig.size()); ++i)
        x_pre[i] = sig[i];

    for (int n = M - 1; n < static_cast<int>(sig.size()); ++n) {
        int idx = 0;
        x_pre[n] = std::complex<double>(0.0, 0.0);

        for (int p = 1; p <= P; p += 2) {
            for (int m = 0; m < M; ++m) {
                x_pre[n] += coeffs[idx++] *
                            sig[n - m] *
                            std::pow(std::abs(sig[n - m]), p - 1);
            }
        }
    }

    return x_pre;
}



