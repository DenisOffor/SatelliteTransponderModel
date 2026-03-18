#include "dpd.h"
#include <cmath>

// ======================
// Вспомогательная функция
// ======================
double DPD::computeRMS(const std::vector<std::complex<double>>& x)
{
    if (x.empty()) return 1.0;

    double sum = 0.0;
    for (const auto& s : x)
        sum += std::norm(s);

    return std::sqrt(sum / x.size());
}

// ======================
// TRAIN (Indirect Learning)
// ======================
void DPD::train(const std::vector<std::complex<double>> &pa_input,
                const std::vector<std::complex<double>> &pa_output,
                int P, int M)
{
    using namespace Eigen;

    size_t N = pa_output.size();
    size_t K = P * M;
    size_t rows = N - M + 1;

    // 🔴 ОДИН scale для всего
    train_scale = computeRMS(pa_input);

    std::vector<std::complex<double>> x_norm(N);
    std::vector<std::complex<double>> y_norm(N);

    for (size_t i = 0; i < N; ++i) {
        x_norm[i] = pa_input[i]  / train_scale;
        y_norm[i] = pa_output[i] / train_scale;
    }

    MatrixXcd Phi(rows, K);
    VectorXcd x(rows);

    // Формирование матрицы
    for (size_t n = M-1; n < N; ++n) {
        size_t row = n - (M-1);

        x(row) = x_norm[n];  // ЦЕЛЬ: восстановить вход PA

        size_t col = 0;

        for (int m = 0; m < M; ++m) {
            std::complex<double> y = y_norm[n - m];

            for (int p = 1; p <= P; ++p) {
                std::complex<double> basis =
                    std::pow(std::abs(y), 2*(p-1)) * y;

                Phi(row, col++) = basis;
            }
        }
    }

    // 🔧 МНК (с минимальной стабилизацией)
    MatrixXcd PhiH = Phi.adjoint();
    MatrixXcd A = (PhiH * Phi + 1e-8 * MatrixXcd::Identity(K, K)).ldlt().solve(PhiH * x);

    coeffs.resize(K);
    for (size_t k = 0; k < K; ++k)
        coeffs[k] = A(k);
}

// ======================
// APPLY DPD
// ======================
std::vector<std::complex<double>> DPD::applyPreDistortion(
    const std::vector<std::complex<double>> &input_signal,
    int P, int M)
{
    size_t N = input_signal.size();
    std::vector<std::complex<double>> x_pre(N, {0.0, 0.0});

    // ✅ Используем ТОТ ЖЕ scale, что в train
    std::vector<std::complex<double>> input_norm(N);
    for (size_t i = 0; i < N; ++i)
        input_norm[i] = input_signal[i];

    // первые отсчёты без изменений
    for (size_t i = 0; i < M-1 && i < N; ++i)
        x_pre[i] = input_norm[i];

    // основной цикл
    for (size_t n = M-1; n < N; ++n) {

        std::complex<double> sum = 0.0;
        size_t col = 0;

        for (int m = 0; m < M; ++m) {

            std::complex<double> y = input_norm[n - m];

            for (int p = 1; p <= P; ++p) {

                std::complex<double> basis =
                    std::pow(std::abs(y), 2*(p-1)) * y;

                sum += coeffs[col++] * basis;
            }
        }

        x_pre[n] = sum;
    }

    // ✅ ВОССТАНАВЛИВАЕМ масштаб
    //for (size_t i = 0; i < N; ++i)
        //x_pre[i] *= train_scale;

    return x_pre;
}
