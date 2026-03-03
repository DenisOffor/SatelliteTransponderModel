#include "pamodels.h"

PAModels::PAModels() {}

void PAModels::SalehModel(std::vector<std::complex<double>>& sig,
                          std::vector<double>& Coeffs,
                          int& linear_gain_dB, int& IBO_dB)
{
    std::vector<double> amplitude_in(sig.size());
    std::vector<double> phase_in(sig.size());
    std::vector<double> amplitude_out(sig.size());
    std::vector<double> phase_out(sig.size());

    double gain_linear = qPow(10, linear_gain_dB / 10.0);

    double A_sat = 1.0 / std::sqrt(Coeffs[1]);
    apply_IBO(sig, IBO_dB, A_sat);

    for(int i = 0; i < sig.size(); ++i) {
        amplitude_in[i] = std::abs(sig[i]);
        phase_in[i] = std::arg(sig[i]);

        // Модель Saleh для АМ/АМ и АМ/ФМ характеристик
        amplitude_out[i] = (Coeffs[0] * amplitude_in[i]) /
                           (1 + Coeffs[1] * amplitude_in[i] * amplitude_in[i]) *
                           gain_linear;

        phase_out[i] = (Coeffs[2] * amplitude_in[i] * amplitude_in[i]) /
                       (1 + Coeffs[3] * amplitude_in[i] * amplitude_in[i]);

        // Применяем искажения
        sig[i] = std::polar(amplitude_out[i], phase_in[i] + phase_out[i]);
    }
}

void PAModels::RappModel(std::vector<std::complex<double>>& sig, std::vector<double>& Coeffs, int& linear_gain_dB, int& IBO_dB)
{
    std::vector<double> amplitude_in(sig.size());
    std::vector<double> phase_in(sig.size());
    std::vector<double> amplitude_out(sig.size());
    std::vector<double> phase_out(sig.size());

    double gain_linear = qPow(10, linear_gain_dB / 10.0);

    double A_sat = Coeffs[0];
    apply_IBO(sig, IBO_dB, A_sat);
    for(int i = 0; i < sig.size(); ++i) {
        amplitude_in[i] = std::abs(sig[i]);
        phase_in[i] = std::arg(sig[i]);

        amplitude_out[i] = gain_linear * amplitude_in[i] /
                           qPow(1 + qPow(amplitude_in[i] / Coeffs[0], 2 * Coeffs[1]), 1 / (2 * Coeffs[1]));
        phase_out[i] = 0;

        // Применяем искажения
        sig[i] = std::polar(amplitude_out[i], phase_in[i] + phase_out[i]);
    }
}

void PAModels::GhorbaniModel(std::vector<std::complex<double>>& sig, std::vector<double>& Coeffs, int& linear_gain_dB, int& IBO_dB)
{
    std::vector<double> amplitude_in(sig.size());
    std::vector<double> phase_in(sig.size());
    std::vector<double> amplitude_out(sig.size());
    std::vector<double> phase_out(sig.size());

    double gain_linear = qPow(10, linear_gain_dB / 10.0);
    double A_sat = find_Asat_Ghorbani(Coeffs, gain_linear);
    apply_IBO(sig, IBO_dB, A_sat);
    double A_rms = 0.0;

    for(int i = 0; i < sig.size(); ++i) {
        amplitude_in[i] = std::abs(sig[i]);
        phase_in[i] = std::arg(sig[i]);

        amplitude_out[i] = gain_linear * Coeffs[0] * amplitude_in[i]
                           / (1 + Coeffs[1] * qPow(amplitude_in[i], 2) + Coeffs[2] *
                                                                             qPow(amplitude_in[i], 4));;
        phase_out[i] = Coeffs[3] * qPow(amplitude_in[i], 2) / (1 + Coeffs[4]
                                                                       * qPow(amplitude_in[i], 2) + Coeffs[5]
                                                                     * qPow(amplitude_in[i], 4));

        // Применяем искажения
        sig[i] = std::polar(amplitude_out[i], phase_in[i] + phase_out[i]);
    }
}

void PAModels::WienerModel(std::vector<std::complex<double>>& sig, QString Static_model,
                           std::vector<double>& Coeffs, std::vector<double>& FIR_Coeffs,
                           int& linear_gain_dB, int& IBO_dB)
{
    // Применяем статическую нелинейность
    if(Static_model == "Saleh")
        SalehModel(sig, Coeffs, linear_gain_dB, IBO_dB);
    else if (Static_model == "Rapp")
        RappModel(sig, Coeffs, linear_gain_dB, IBO_dB);
    else if (Static_model == "Ghorbani")
        GhorbaniModel(sig, Coeffs, linear_gain_dB, IBO_dB);

    ApplyFIRWithMemory(sig, FIR_Coeffs, 5);  // 3 - лучше сделать параметром
}

void PAModels::apply_IBO(std::vector<std::complex<double>>& tx,
                         double IBO_dB,
                         double A_sat)
{
    if (tx.empty())
        return;

    double sum = 0.0;
    for (const auto& s : tx)
        sum += std::norm(s);

    double rms = std::sqrt(sum / tx.size());

    if (rms == 0.0)
        return;

    // Входная мощность относительно насыщения
    double P_sat = A_sat * A_sat;

    double P_in_target = P_sat / std::pow(10.0, IBO_dB / 10.0);

    double A_rms_target = std::sqrt(P_in_target);

    double scale = A_rms_target / rms;

    for (auto& s : tx)
        s *= scale;
}

double PAModels::find_Asat_Ghorbani(const std::vector<double>& c,
                                    double gain_linear)
{
    double Amax = 10.0;      // верхняя граница поиска
    double step = 0.0005;    // шаг

    double A_sat = 0.0;
    double max_out = 0.0;

    for (double Ain = 0.0; Ain < Amax; Ain += step)
    {
        double Aout = c[0] * Ain /
                      (1.0 + c[1]*Ain*Ain + c[2]*pow(Ain,4));

        if (Aout > max_out)
        {
            max_out = Aout;
            A_sat = Ain;
        }
    }

    return A_sat;
}

void PAModels::ApplyFIRWithMemory(std::vector<std::complex<double>>& signal, const std::vector<double>& FIR_Coefs, int numTaps)
{
    const size_t N = signal.size();

    if (N == 0 || FIR_Coefs.size() < 2 || numTaps <= 0) {
        return;  // Проверка входных данных
    }

    double C = FIR_Coefs[0];
    double alpha = FIR_Coefs[1];

    // Формируем коэффициенты FIR
    std::vector<std::complex<double>> h(numTaps);
    for (int m = 0; m < numTaps; ++m)
    {
        h[m] = C * std::pow(alpha, m);  // Вещественные коэффициенты
    }

    // Нормировка
    double sum = 0.0;
    for (int m = 0; m < numTaps; ++m)
        sum += std::abs(h[m]);  // или real(h[m]) т.к. они вещественные

    if (sum > 0) {
        for (int m = 0; m < numTaps; ++m)
            h[m] /= sum;
    }

    // Буфер для результата
    std::vector<std::complex<double>> output(N, {0.0, 0.0});

    // FIR фильтрация
    for (size_t n = 0; n < N; ++n)
    {
        for (int m = 0; m < numTaps; ++m)
        {
            if (n >= static_cast<size_t>(m))
            {
                output[n] += h[m] * signal[n - m];
            }
        }
    }

    // Копируем результат обратно
    signal = std::move(output);
}

