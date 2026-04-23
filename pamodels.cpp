#include "pamodels.h"

PAModels::PAModels() {}

double computeAvgPower(const std::vector<std::complex<double>>& x)
{
    if (x.empty()) return 0.0;

    double sum = 0.0;
    for (const auto& s : x)
        sum += std::norm(s);

    return sum / x.size();
}

double computeMeasuredIBO_dB(const std::vector<std::complex<double>>& x, double A_sat)
{
    double Pin_avg = computeAvgPower(x);
    double Psat_in = A_sat * A_sat;

    if (Pin_avg <= 0.0) return 0.0;

    return 10.0 * std::log10(Psat_in / Pin_avg);
}

void PAModels::SalehModel(std::vector<std::complex<double>>& sig,
                          std::vector<double>& Coeffs,
                          int& linear_gain_dB, int& IBO_dB)
{
    std::vector<double> amplitude_in(sig.size());
    std::vector<double> phase_in(sig.size());
    std::vector<double> amplitude_out(sig.size());
    std::vector<double> phase_out(sig.size());

    double gain_linear = qPow(10, linear_gain_dB / 20.0);

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

    double gain_linear = qPow(10, linear_gain_dB / 20.0);

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

    double gain_linear = qPow(10, linear_gain_dB / 20.0);

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
    ApplyFIRWithMemory(sig, FIR_Coeffs[0], 9);  // 3 - лучше сделать параметром
    // Применяем статическую нелинейность
    if(Static_model == "Saleh")
        SalehModel(sig, Coeffs, linear_gain_dB, IBO_dB);
    else if (Static_model == "Rapp")
        RappModel(sig, Coeffs, linear_gain_dB, IBO_dB);
    else if (Static_model == "Ghorbani")
        GhorbaniModel(sig, Coeffs, linear_gain_dB, IBO_dB);
}

void PAModels::HammersteinModel(std::vector<std::complex<double>> &sig, QString Static_model, std::vector<double> &Coeffs, std::vector<double> &FIR_Coeffs, int &linear_gain_dB, int &IBO_dB)
{
    // Применяем статическую нелинейность
    if(Static_model == "Saleh")
        SalehModel(sig, Coeffs, linear_gain_dB, IBO_dB);
    else if (Static_model == "Rapp")
        RappModel(sig, Coeffs, linear_gain_dB, IBO_dB);
    else if (Static_model == "Ghorbani")
        GhorbaniModel(sig, Coeffs, linear_gain_dB, IBO_dB);

    ApplyFIRWithMemory(sig, FIR_Coeffs[0], 9);
}

void PAModels::WHModel(std::vector<std::complex<double> > &sig, QString Static_model, std::vector<double> &Coeffs, std::vector<double> &FIR_Coeffs, int &linear_gain_dB, int &IBO_dB)
{
    ApplyFIRWithMemory(sig, FIR_Coeffs[0], 9);

    // Применяем статическую нелинейность
    if(Static_model == "Saleh")
        SalehModel(sig, Coeffs, linear_gain_dB, IBO_dB);
    else if (Static_model == "Rapp")
        RappModel(sig, Coeffs, linear_gain_dB, IBO_dB);
    else if (Static_model == "Ghorbani")
        GhorbaniModel(sig, Coeffs, linear_gain_dB, IBO_dB);

    ApplyFIRWithMemory(sig, FIR_Coeffs[1], 9);
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

void PAModels::ApplyFIRWithMemory(std::vector<std::complex<double>>& signal, double alpha, int numTaps)
{
    const size_t N = signal.size();

    // Формируем коэффициенты FIR
    std::vector<std::complex<double>> h(numTaps);
    for (int m = 0; m < numTaps; ++m)
    {
        double mag = std::pow(alpha, m);
        double phase = alpha / 3 * m;
        h[m] = std::polar(mag, phase);
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

void PAModels::ScaleToRMS_forPA(std::vector<std::complex<double>>& sig, Source& source)
{
    double A_sat;
    double target_rms;
    if (source.PAModel == "Saleh")
        A_sat = 1.0 / std::sqrt(source.SalehCoeffs[1]);
    else if(source.PAModel == "Rapp")
        A_sat = source.RappCoeffs[0];
    else if(source.PAModel == "Ghorbani")
        A_sat = find_Asat_Ghorbani(source.GhorbaniCoeffs, qPow(10, source.linear_gain_dB / 20.0));
    else if(source.PAModel == "Wiener") {
        if (source.W_StaticNonlinModel == "Saleh")
            A_sat = 1.0 / std::sqrt(source.SalehCoeffs[1]);
        else if(source.W_StaticNonlinModel == "Rapp")
            A_sat = source.RappCoeffs[0];
        else if(source.W_StaticNonlinModel == "Ghorbani")
            A_sat = find_Asat_Ghorbani(source.GhorbaniCoeffs, qPow(10, source.linear_gain_dB / 20.0));
    }
    else if(source.PAModel == "Hammerstein") {
        if (source.H_StaticNonlinModel == "Saleh")
            A_sat = 1.0 / std::sqrt(source.SalehCoeffs[1]);
        else if(source.H_StaticNonlinModel == "Rapp")
            A_sat = source.RappCoeffs[0];
        else if(source.H_StaticNonlinModel == "Ghorbani")
            A_sat = find_Asat_Ghorbani(source.GhorbaniCoeffs, qPow(10, source.linear_gain_dB / 20.0));
    }
    else if(source.PAModel == "Wiener-Hammerstein") {
        if (source.WH_StaticNonlinModel == "Saleh")
            A_sat = 1.0 / std::sqrt(source.SalehCoeffs[1]);
        else if(source.WH_StaticNonlinModel == "Rapp")
            A_sat = source.RappCoeffs[0];
        else if(source.WH_StaticNonlinModel == "Ghorbani")
            A_sat = find_Asat_Ghorbani(source.GhorbaniCoeffs, qPow(10, source.linear_gain_dB / 20.0));
    }
    source.PA_sat = A_sat * A_sat * qPow(10, source.linear_gain_dB / 10.0);
    target_rms = A_sat / std::pow(10.0, source.IBO_dB / 20.0);

    scaleToRMS(sig, target_rms);
}

void PAModels::scaleToRMS(std::vector<std::complex<double>>& x, double target_rms)
{
    double sum = 0.0;
    for (const auto& s : x) sum += std::norm(s);
    double rms = std::sqrt(sum / x.size());
    if (rms == 0.0) return;

    double k = target_rms / rms;
    for (auto& s : x) s *= k;
}

