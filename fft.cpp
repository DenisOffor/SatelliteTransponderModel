#include "fft.h"
#include <QDebug>

FFT::FFT(int n) : m_n(n)
{
    m_buffer.resize(n);

    m_planForward = fftw_plan_dft_1d(
        n,
        reinterpret_cast<fftw_complex*>(m_buffer.data()),
        reinterpret_cast<fftw_complex*>(m_buffer.data()),
        FFTW_FORWARD,
        FFTW_MEASURE
        );

    m_planBackward = fftw_plan_dft_1d(
        n,
        reinterpret_cast<fftw_complex*>(m_buffer.data()),
        reinterpret_cast<fftw_complex*>(m_buffer.data()),
        FFTW_BACKWARD,
        FFTW_MEASURE
        );
}

FFT::~FFT()
{
    if(m_plansCreated) {
        fftw_destroy_plan(m_planForward);
        fftw_destroy_plan(m_planBackward);
    }
}

// std::vector<std::complex<double>> FFT::fft(const std::vector<std::complex<double>>& x)
// {
//     std::vector<std::complex<double>> result(m_n);

//     if(x.size() != m_n || !m_plansCreated) {
//         qWarning() << "FFTWWrapper: размер не совпадает или планы не созданы";
//         return result;
//     }

//     // Копируем входные данные
//     for(int i = 0; i < m_n; ++i) {
//         m_in[i][0] = x[i].real();
//         m_in[i][1] = x[i].imag();
//     }

//     // Выполняем FFT
//     fftw_execute(m_planForward);

//     // Копируем результат
//     for(int i = 0; i < m_n; ++i) {
//         result[i] = std::complex<double>(m_out[i][0], m_out[i][1]);
//     }

//     return result;
// }

std::vector<std::complex<double>> FFT::fftInPlace(std::vector<std::complex<double>> x)
{
    fftw_execute_dft(
        m_planForward,
        reinterpret_cast<fftw_complex*>(x.data()),
        reinterpret_cast<fftw_complex*>(x.data())
        );
    return x;
}

std::vector<std::complex<double>> FFT::ifftInPlace(std::vector<std::complex<double>> x)
{
    fftw_execute_dft(
        m_planBackward,
        reinterpret_cast<fftw_complex*>(x.data()),
        reinterpret_cast<fftw_complex*>(x.data())
        );
    return x;
}

// std::vector<std::complex<double>> FFT::ifft(const std::vector<std::complex<double>>& X)
// {
//     std::vector<std::complex<double>> result(m_n);

//     if(X.size() != m_n || !m_plansCreated) {
//         qWarning() << "FFTWWrapper: размер не совпадает или планы не созданы";
//         return result;
//     }

//     // Копируем входные данные
//     for(int i = 0; i < m_n; ++i) {
//         m_in[i][0] = X[i].real();
//         m_in[i][1] = X[i].imag();
//     }

//     // Выполняем обратное FFT
//     fftw_execute(m_planBackward);

//     // Копируем результат и делим на N
//     for(int i = 0; i < m_n; ++i) {
//         result[i] = std::complex<double>(m_out[i][0] / m_n, m_out[i][1] / m_n);
//     }

//     return result;
// }
