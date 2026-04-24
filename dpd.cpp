#include "dpd.h"
#include <cmath>

void DPD::train(const std::vector<std::complex<double>>& pa_input,
                const std::vector<std::complex<double>>& pa_output,
                const Source& source)
{
    int P = 0;
    int M = 0;

    if (source.PredistorterType == "MP") {
        P = source.MP_P;
        M = source.MP_M;
    }
    else if (source.PredistorterType == "GMP") {
        P = source.GMP_P;
        M = source.GMP_M;
    }

    std::vector<std::complex<double>> pa_output_norm;
    if (source.NormalizationType == "Peak normalization") {
        Gpeak = computePeakGain(pa_input, pa_output);
        pa_output_norm = normalizeByGain(pa_output, Gpeak);
    }
    else if (source.NormalizationType == "RMS normalization") {
        Gpeak = computePeakGain(pa_input, pa_output);
        Grms = computeGrms(pa_input, pa_output);
        pa_output_norm = normalizeByGain(pa_output, Grms);
    }
    else if(source.NormalizationType == "Complex linear gain") {
        Gpeak = computeLinearGainLSAbs(pa_input, pa_output);
        pa_output_norm = normalizeByGain(pa_output, Gpeak);
    }

    QElapsedTimer timer;
    timer.start();

    VectorXcd a;
    int num_coeffs = 0;

    if (source.PredistorterType == "MP") {
        a = DPDsolve_least_squares_fixed_linear(
            make_MP_mat(pa_output_norm, P, M, source.Enable_even_P),
            make_goal(pa_input, source)
            );
    }
    else if (source.PredistorterType == "GMP") {
        a = DPDsolve_least_squares_fixed_linear(
            make_GMP_mat(pa_output_norm, P, M,
                         source.GMP_L_lag, source.GMP_L_lead, source.Enable_even_P),
            make_goal(pa_input, source)
            );
    }

    qDebug() << "LS train time:" << timer.elapsed() << "ms";

    coeffs.resize(a.size());
    for (int i = 0; i < a.size(); ++i)
        coeffs[i] = a(i);
}

double DPD::computeLinearGainLSAbs(const std::vector<std::complex<double>>& pa_input, const std::vector<std::complex<double>>& pa_output)
{
    return std::abs(computeLinearGainLS(pa_input, pa_output));
}

std::complex<double> DPD::computeLinearGainLS(
    const std::vector<std::complex<double>>& pa_input,
    const std::vector<std::complex<double>>& pa_output)
{
    if (pa_input.empty() || pa_output.empty() || pa_input.size() != pa_output.size())
        return {1.0, 0.0};

    std::complex<double> num(0.0, 0.0);
    double den = 0.0;

    for (size_t i = 0; i < pa_input.size(); ++i)
    {
        num += std::conj(pa_input[i]) * pa_output[i];
        den += std::norm(pa_input[i]);
    }

    if (den <= 1e-15)
        return {1.0, 0.0};

    return num / den;
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

MatrixXcd DPD::make_MP_mat(const std::vector<std::complex<double>>& x, const int P, const int M, bool Enable_even_P)
{
    // y = Phi(x...)*a
    // исключили первые M-1 отсчетов, чтобы не уходить в x[-1] и тд, то есть отрицательные отсчеты
    int N = x.size();
    int rows = N - M + 1;
    int cols;
    int p_iter;
    if(Enable_even_P) { cols = M * P; p_iter = 1;}
    else { cols = M * (P + 1) / 2; p_iter = 2; }

    MatrixXcd Phi_LS(rows, cols);
    qDebug() << "MP mat: " << rows << "x" << cols;

    for(int n = M - 1; n < N; ++n) {
        size_t cur_row = n - M + 1;
        size_t col_idx = 0;

        for(int p = 1; p <= P; p += p_iter)
            for(int m = 0; m < M; ++m)
                Phi_LS(cur_row, col_idx++) = x[n - m] * std::pow(std::abs(x[n - m]), p - 1);
    }

    return Phi_LS;
}

MatrixXcd DPD::make_GMP_mat(const std::vector<std::complex<double>>& x,
                            const int P, const int M,
                            const int L_lag, const int L_lead, bool Enable_even_P)
{
    int N = x.size();
    int rows = N - M - L_lag - L_lead + 1;
    int cols;
    int p_iter;
    if(Enable_even_P) { cols = M * (1 + L_lag + L_lead) * P; p_iter = 1; }
    else { cols = M * (1 + L_lag + L_lead) * ((P + 1) / 2); p_iter = 2; }
    MatrixXcd Phi_LS(rows, cols);

    qDebug() << "GMP mat:" << rows << "x" << cols;

    for (int n = M - 1 + L_lag; n < N - L_lead; ++n) {
        int cur_row = n - (M - 1 + L_lag);
        int col_idx = 0;

        // aligned
        for (int p = 1; p <= P; p += p_iter)
            for (int m = 0; m < M; ++m)
                Phi_LS(cur_row, col_idx++) =
                    x[n - m] * std::pow(std::abs(x[n - m]), p - 1);

        // lagging
        for (int p = 1; p <= P; p += p_iter)
            for (int m = 0; m < M; ++m)
                for (int l = 1; l <= L_lag; ++l)
                    Phi_LS(cur_row, col_idx++) =
                        x[n - m] * std::pow(std::abs(x[n - m - l]), p - 1);

        // leading
        for (int p = 1; p <= P; p += p_iter)
            for (int m = 0; m < M; ++m)
                for (int l = 1; l <= L_lead; ++l)
                    Phi_LS(cur_row, col_idx++) =
                        x[n - m] * std::pow(std::abs(x[n - m + l]), p - 1);
    }

    return Phi_LS;
}

VectorXcd DPD::make_goal(const std::vector<std::complex<double>>& y, const Source& source) {
    int N = y.size();
    int n_start = 0;
    int n_end = 0;
    if(source.PredistorterType == "MP") {
        n_start = source.MP_M - 1;
        n_end = N;
    }
    else if (source.PredistorterType == "GMP") {
        n_start = source.GMP_M - 1 + source.GMP_L_lag;
        n_end = N - source.GMP_L_lead;
    }
    VectorXcd y_LS(n_end - n_start);

    size_t num_el = 0;
    for(int n = n_start; n < n_end; ++n)
        y_LS(num_el++) = y[n];

    return y_LS;
}

VectorXcd DPD::DPDsolve_least_squares(const MatrixXcd& Phi, const VectorXcd& goal) {
    // Нахождение решения МНК через нормальное уравнение
    //return (Phi.adjoint() * Phi).ldlt().solve(Phi.adjoint() * goal);
    return Phi.colPivHouseholderQr().solve(goal);
}

VectorXcd DPD::DPDsolve_least_squares_fixed_linear(const MatrixXcd &A, const VectorXcd &b)
{
    if (A.cols() == 0 || A.rows() == 0)
        return VectorXcd();

    // Первый коэффициент фиксируем
    const std::complex<double> a0_fixed(1.0, 0.0);

    // Если в модели только один коэффициент
    if (A.cols() == 1) {
        VectorXcd a(1);
        a(0) = a0_fixed;
        return a;
    }

    // Первый столбец = главный линейный член
    VectorXcd phi0 = A.col(0);

    // Новый целевой вектор:
    // b = phi0 * 1 + Arest * arest
    VectorXcd b_rest = b - phi0 * a0_fixed;

    // Матрица без первого столбца
    MatrixXcd A_rest = A.rightCols(A.cols() - 1);

    // Решаем для остальных коэффициентов
    VectorXcd a_rest =
        A_rest.colPivHouseholderQr().solve(b_rest); //(A_rest.adjoint() * A_rest).ldlt().solve(A_rest.adjoint() * b_rest);

    // Склеиваем полный вектор коэффициентов
    VectorXcd a(A.cols());
    a(0) = a0_fixed;
    a.tail(A.cols() - 1) = a_rest;

    return a;
}

std::vector<std::complex<double>> DPD::applyMP(
    const std::vector<std::complex<double>>& sig, const Source& source)
{
    int P = source.MP_P;
    int M = source.MP_M;
    std::vector<std::complex<double>> x_pre(sig.size(), std::complex<double>(0.0, 0.0));

    double norm_coef = 1.0;
    if(source.NormalizationType == "RMS normalization")
    //    norm_coef = Gpeak/Grms;

    for (int i = 0; i < M - 1 && i < static_cast<int>(sig.size()); ++i)
        x_pre[i] = sig[i] * norm_coef;

    int p_iter = source.Enable_even_P ? 1 : 2;

    for (int n = M - 1; n < static_cast<int>(sig.size()); ++n) {
        int idx = 0;
        x_pre[n] = std::complex<double>(0.0, 0.0);

        for (int p = 1; p <= P; p += p_iter) {
            for (int m = 0; m < M; ++m) {
                x_pre[n] += coeffs[idx++] *
                            sig[n - m] * norm_coef *
                            std::pow(std::abs(sig[n - m] * norm_coef), p - 1);
            }
        }
    }

    return x_pre;
}

std::vector<std::complex<double>> DPD::applyGMP(
    const std::vector<std::complex<double>>& sig, const Source& source)
{
    int P = source.GMP_P;
    int M = source.GMP_M;
    int L_lag = source.GMP_L_lag;
    int L_lead = source.GMP_L_lead;
    std::vector<std::complex<double>> x_pre(sig.size(), std::complex<double>(0.0, 0.0));

    double norm_coef = 1.0;
    if(source.NormalizationType == "RMS normalization")
    //    norm_coef = Gpeak/Grms;

    for (int i = 0; i < M - 1 + L_lag && i < static_cast<int>(sig.size()); ++i)
        x_pre[i] = sig[i] * norm_coef;

    int p_iter = source.Enable_even_P ? 1 : 2;

    for (int n = M - 1 + L_lag; n < static_cast<int>(sig.size()) - L_lead; ++n) {
        int idx = 0;
        x_pre[n] = std::complex<double>(0.0, 0.0);

        // 1) aligned
        for (int p = 1; p <= P; p += p_iter) {
            for (int m = 0; m < M; ++m) {
                x_pre[n] += coeffs[idx++] *
                            sig[n - m] * norm_coef *
                            std::pow(std::abs(sig[n - m] * norm_coef), p - 1);
            }
        }

        // 2) lagging
        for (int p = 1; p <= P; p += p_iter) {
            for (int m = 0; m < M; ++m) {
                for (int l = 1; l <= L_lag; ++l) {
                    x_pre[n] += coeffs[idx++] *
                                sig[n - m] * norm_coef *
                                std::pow(std::abs(sig[n - m - l] * norm_coef), p - 1);
                }
            }
        }

        // 3) leading
        for (int p = 1; p <= P; p += p_iter) {
            for (int m = 0; m < M; ++m) {
                for (int l = 1; l <= L_lead; ++l) {
                    x_pre[n] += coeffs[idx++] *
                                sig[n - m] * norm_coef *
                                std::pow(std::abs(sig[n - m + l] * norm_coef), p - 1);
                }
            }
        }
    }

    for (int i = std::max(0, static_cast<int>(sig.size()) - L_lead);
         i < static_cast<int>(sig.size()); ++i)
        x_pre[i] = sig[i] * norm_coef;

    return x_pre;
}

void DPD::softClip(std::vector<std::complex<double>>& x, double A_lim)
{
    for (auto& s : x) {
        double a = std::abs(s);
        if (a > 1e-12) {
            double gain = std::tanh(a / A_lim) / (a / A_lim);
            s *= gain;
        }
    }
}

double DPD::getAmplitudeLimitFromTargetPAPR(const std::vector<std::complex<double>>& x,
                                       double paprTarget_dB)
{
    double rms = computeRMS(x);
    return rms * std::pow(10.0, paprTarget_dB / 20.0);
}

