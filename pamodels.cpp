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

    apply_IBO(sig, IBO_dB);
    double sum = 0.0;

    for (const auto& sample : sig)
    {
        sum += std::norm(sample);   // |x|^2 = real^2 + imag^2
    }

    double A_rms = std::sqrt(sum / sig.size());

    qDebug() << "A_rms =" << A_rms;
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

    apply_IBO(sig, IBO_dB);
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

    apply_IBO(sig, IBO_dB);
    double sum = 0.0;

    for (const auto& sample : sig)
    {
        sum += std::norm(sample);   // |x|^2 = real^2 + imag^2
    }

    double A_rms = std::sqrt(sum / sig.size());

    qDebug() << "A_rms =" << A_rms;

    for(int i = 0; i < sig.size(); ++i) {
        amplitude_in[i] = std::abs(sig[i]);
        phase_in[i] = std::arg(sig[i]);

        amplitude_out[i] = gain_linear * Coeffs[0] * amplitude_in[i]
                           / (1 + Coeffs[1] * qPow(amplitude_in[i], 2) + Coeffs[2] *
                            qPow(amplitude_in[i], 4));;
        phase_out[i] = Coeffs[3] * qPow(amplitude_in[i], 2) / (1 + Coeffs[4]
                            * qPow(amplitude_in[i], 2) + Coeffs[5]
                            * qPow(amplitude_in[i], 4)) * 180 / 3.14;

        sig[i] = std::polar(amplitude_out[i], phase_in[i] + phase_out[i]);
    }
}

void PAModels::apply_IBO(QVector<std::complex<double>>& tx, double IBO_dB)
{
    if (tx.isEmpty()) {
        return;
    }

    // Вычисляем RMS
    double sum_squares = 0.0;
    for (const auto& sample : tx) {
        sum_squares += std::norm(sample);
    }
    double tx_rms = std::sqrt(sum_squares / tx.size());

    // Защита от деления на ноль
    if (tx_rms == 0.0) {
        return;
    }

    // Вычисляем целевое значение RMS с учетом IBO
    double A_rms_target = std::pow(10.0, -IBO_dB / 20.0);

    // Применяем нормировку и IBO за один проход
    double scale = A_rms_target / tx_rms;  // объединенный коэффициент масштабирования

    for (int i = 0; i < tx.size(); ++i) {
        tx[i] = tx[i] * scale;
    }
}
