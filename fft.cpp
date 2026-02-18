#include "fft.h"
#include <QDebug>

FFT::FFT(int n) : m_n(n), m_plansCreated(false)
{
    m_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    m_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);

    if(m_in && m_out) {
        m_planForward = fftw_plan_dft_1d(n, m_in, m_out, FFTW_FORWARD, FFTW_MEASURE);
        m_planBackward = fftw_plan_dft_1d(n, m_in, m_out, FFTW_BACKWARD, FFTW_MEASURE);
        m_plansCreated = true;
    }
}

FFT::~FFT()
{
    if(m_plansCreated) {
        fftw_destroy_plan(m_planForward);
        fftw_destroy_plan(m_planBackward);
    }
    fftw_free(m_in);
    fftw_free(m_out);
}

QVector<std::complex<double>> FFT::fft(const QVector<std::complex<double>>& x)
{
    QVector<std::complex<double>> result(m_n);

    if(x.size() != m_n || !m_plansCreated) {
        qWarning() << "FFTWWrapper: размер не совпадает или планы не созданы";
        return result;
    }

    // Копируем входные данные
    for(int i = 0; i < m_n; ++i) {
        m_in[i][0] = x[i].real();
        m_in[i][1] = x[i].imag();
    }

    // Выполняем FFT
    fftw_execute(m_planForward);

    // Копируем результат
    for(int i = 0; i < m_n; ++i) {
        result[i] = std::complex<double>(m_out[i][0], m_out[i][1]);
    }

    return result;
}

QVector<std::complex<double>> FFT::ifft(const QVector<std::complex<double>>& X)
{
    QVector<std::complex<double>> result(m_n);

    if(X.size() != m_n || !m_plansCreated) {
        qWarning() << "FFTWWrapper: размер не совпадает или планы не созданы";
        return result;
    }

    // Копируем входные данные
    for(int i = 0; i < m_n; ++i) {
        m_in[i][0] = X[i].real();
        m_in[i][1] = X[i].imag();
    }

    // Выполняем обратное FFT
    fftw_execute(m_planBackward);

    // Копируем результат и делим на N
    for(int i = 0; i < m_n; ++i) {
        result[i] = std::complex<double>(m_out[i][0] / m_n, m_out[i][1] / m_n);
    }

    return result;
}
