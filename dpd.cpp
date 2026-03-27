#include "dpd.h"
#include <cmath>

void DPD::train(const std::vector<std::complex<double>> &pa_input,
                const std::vector<std::complex<double> > &pa_output, const Source& source)
{
    int P = source.MP_P;
    int M = source.MP_M;

    std::vector<std::complex<double>> pa_output_norm;
    if(source.NormalizationType == "Peak normalization") {
        Gpeak = computePeakGain(pa_input, pa_output);
        pa_output_norm = normalizeByGain(pa_output, Gpeak);
    }
    else if(source.NormalizationType == "RMS normalization") {
        Gpeak = computePeakGain(pa_input, pa_output);
        Grms = computeGrms(pa_input, pa_output);
        pa_output_norm = normalizeByGain(pa_output, Grms);
    }
    QElapsedTimer timer;
    timer.start();

    VectorXcd a = DPDsolve_least_squares(
        make_MP_mat(pa_output_norm, P, M),
        make_goal(pa_input, M)
        );
    qDebug() << "LS train time:" << timer.elapsed() << "ms";

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

std::vector<std::complex<double>> DPD::normalizeByGain(const std::vector<std::complex<double>>& pa_output,
    double G)
{
    std::vector<std::complex<double>> y_norm(pa_output.size());

    if (G == 0.0)
        return y_norm;

    for (size_t i = 0; i < pa_output.size(); ++i)
        y_norm[i] = pa_output[i] / G;

    return y_norm;
}

double DPD::computeRMS(const std::vector<std::complex<double>>& sig)
{
    if (sig.empty()) return 0.0;

    double sum = 0.0;
    for (const auto& s : sig)
        sum += std::norm(s);

    return std::sqrt(sum / sig.size());
}

double DPD::computeGrms(const std::vector<std::complex<double>>& pa_input,
                   const std::vector<std::complex<double>>& pa_output)
{
    double rms_in = computeRMS(pa_input);
    double rms_out = computeRMS(pa_output);

    if (rms_in == 0.0)
        return 1.0;

    return rms_out / rms_in;
}

MatrixXcd DPD::make_MP_mat(const std::vector<std::complex<double>>& x, const int P, const int M)
{
    const int N = static_cast<int>(x.size());
    const int rows = N - M + 1;
    const int cols = M * ((P + 1) / 2);
    const int L = (P + 1) / 2;

    MatrixXcd Phi_LS(rows, cols);

    for (int n = M - 1; n < N; ++n) {
        const int cur_row = n - M + 1;
        int col_idx = 0;

        for (int p_idx = 0; p_idx < L; ++p_idx) {
            for (int m = 0; m < M; ++m) {
                const std::complex<double> s = x[n - m];
                const double a2 = std::norm(s);

                std::complex<double> term = s;
                for (int k = 0; k < p_idx; ++k)
                    term *= a2;

                Phi_LS(cur_row, col_idx++) = term;
            }
        }
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
    const std::vector<std::complex<double>>& sig, const Source& source)
{
    int P = source.MP_P;
    int M = source.MP_M;
    std::vector<std::complex<double>> x_pre(sig.size(), std::complex<double>(0.0, 0.0));

    double norm_coef = 1.0;
    if(source.NormalizationType == "RMS normalization")
        norm_coef = Gpeak/Grms;

    for (int i = 0; i < M - 1 && i < static_cast<int>(sig.size()); ++i)
        x_pre[i] = sig[i] * norm_coef;

    for (int n = M - 1; n < static_cast<int>(sig.size()); ++n) {
        int idx = 0;
        x_pre[n] = std::complex<double>(0.0, 0.0);

        for (int p = 1; p <= P; p += 2) {
            for (int m = 0; m < M; ++m) {
                x_pre[n] += coeffs[idx++] *
                            sig[n - m] * norm_coef *
                            std::pow(std::abs(sig[n - m] * norm_coef), p - 1);
            }
        }
    }

    return x_pre;
}



