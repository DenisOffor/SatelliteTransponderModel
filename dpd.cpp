#include "dpd.h"

// Нормализация сигнала, сохраняем scale внутри класса
std::vector<std::complex<double>> DPD::normalizeSignal(const std::vector<std::complex<double>> &sig)
{
    size_t N = sig.size();
    std::vector<std::complex<double>> sig_norm(N);

    double max_val = 0.0;
    for (size_t n = 0; n < N; ++n)
        if (std::abs(sig[n]) > max_val) max_val = std::abs(sig[n]);

    scale = (max_val > 0.0) ? max_val : 1.0;

    for (size_t n = 0; n < N; ++n)
        sig_norm[n] = sig[n] / scale;

    return sig_norm;
}

// Обучение DPD
void DPD::train(const std::vector<std::complex<double>> &pa_input,
                const std::vector<std::complex<double>> &pa_output,
                int P, int M)
{
    using namespace Eigen;

    auto pa_input_norm  = normalizeSignal(pa_input);
    auto pa_output_norm = normalizeSignal(pa_output);

    size_t N = pa_output_norm.size();
    size_t K = P * M;  // количество колонок Φ
    size_t rows = N - M + 1; // строки Φ, начинаем с n=M-1

    MatrixXcd Phi(rows, K);
    VectorXcd x(rows);

    // Формируем Φ и x
    for (size_t n = M-1; n < N; ++n) {
        size_t row = n - (M-1);
        x(row) = pa_input_norm[n];  // целевой сигнал (вход PA)

        size_t col = 0;
        for (int m = 0; m < M; ++m) {
            std::complex<double> y = pa_output_norm[n - m];
            std::complex<double> y_pow = 1.0;
            for (int p = 1; p <= P; ++p) {
                y_pow = std::pow(std::abs(y), 2*(p-1)) * y;
                Phi(row, col++) = y_pow;
            }
        }
    }

    // МНК: a = (Φᴴ Φ)^-1 Φᴴ x
    MatrixXcd PhiH = Phi.adjoint();
    MatrixXcd A = (PhiH * Phi).inverse() * (PhiH * x);

    coeffs.resize(K);
    for (size_t k = 0; k < K; ++k)
        coeffs[k] = A(k);
}

// Применение предыскажения
std::vector<std::complex<double>> DPD::applyPreDistortion(
    const std::vector<std::complex<double>> &input_signal,
    int P, int M)
{
    size_t N = input_signal.size();
    size_t K = P * M;
    std::vector<std::complex<double>> x_pre(N, std::complex<double>(0,0));

    // Сначала нормируем входной сигнал тем же scale, что был на обучении
    std::vector<std::complex<double>> input_norm(N);
    for (size_t n = 0; n < N; ++n)
        input_norm[n] = input_signal[n] / scale;

    // Первые M-1 отсчётов копируем
    for (size_t i = 0; i < M-1; ++i)
        x_pre[i] = input_norm[i];

    // Основной цикл с полиномом памяти
    for (size_t n = M-1; n < N; ++n) {
        std::complex<double> sum = 0.0;
        size_t col = 0;
        for (int m = 0; m < M; ++m) {
            std::complex<double> y = input_norm[n - m];
            std::complex<double> y_pow = 1.0;
            for (int p = 1; p <= P; ++p) {
                y_pow = std::pow(std::abs(y), 2*(p-1)) * y;
                sum += coeffs[col++] * y_pow;
            }
        }
        x_pre[n] = sum;
    }

    // Восстанавливаем исходный уровень сигнала
    for (size_t n = 0; n < N; ++n)
        x_pre[n] *= scale;

    return x_pre;
}
