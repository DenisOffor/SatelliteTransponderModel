#include "pamodels.h"

PAModels::PAModels() {}

void PAModels::SalehModel(QVector<std::complex<double>>& sig,
                          QVector<double>& Coeffs,
                          int& linear_gain_dB, int& IBO_dB)
{
    QVector<double> amplitude_in(sig.size());
    QVector<double> phase_in(sig.size());
    QVector<double> amplitude_out(sig.size());
    QVector<double> phase_out(sig.size());

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

void PAModels::RappModel(QVector<std::complex<double>>& sig, QVector<double>& Coeffs, int& linear_gain_dB, int& IBO_dB)
{
    QVector<double> amplitude_in(sig.size());
    QVector<double> phase_in(sig.size());
    QVector<double> amplitude_out(sig.size());
    QVector<double> phase_out(sig.size());

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

void PAModels::GhorbaniModel(QVector<std::complex<double>>& sig, QVector<double>& Coeffs, int& linear_gain_dB, int& IBO_dB)
{
    QVector<double> amplitude_in(sig.size());
    QVector<double> phase_in(sig.size());
    QVector<double> amplitude_out(sig.size());
    QVector<double> phase_out(sig.size());

    double gain_linear = qPow(10, linear_gain_dB / 10.0);
    double A_sat = find_Asat_Ghorbani(Coeffs, gain_linear);
    apply_IBO(sig, IBO_dB, A_sat);
    double A_rms = 0.0;

    if (!sig.isEmpty())
    {
        double sum = 0.0;

        for (const auto& sample : sig)
            sum += std::norm(sample);

        A_rms = std::sqrt(sum / sig.size());
    }

    qDebug() << "A_rms =" << A_rms;
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

void PAModels::apply_IBO(QVector<std::complex<double>>& tx,
                         double IBO_dB,
                         double A_sat)
{
    if (tx.isEmpty())
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

double PAModels::find_Asat_Ghorbani(const QVector<double>& c,
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
