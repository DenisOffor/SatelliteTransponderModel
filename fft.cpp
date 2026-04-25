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

    m_plansCreated = (m_planForward != nullptr && m_planBackward != nullptr);
}

FFT::~FFT()
{
    if(m_plansCreated) {
        fftw_destroy_plan(m_planForward);
        fftw_destroy_plan(m_planBackward);
    }
}

void FFT::fftInPlace(std::vector<std::complex<double>>& x)
{
    fftw_execute_dft(
        m_planForward,
        reinterpret_cast<fftw_complex*>(x.data()),
        reinterpret_cast<fftw_complex*>(x.data())
        );
}

void FFT::ifftInPlace(std::vector<std::complex<double>>& x)
{
    fftw_execute_dft(
        m_planBackward,
        reinterpret_cast<fftw_complex*>(x.data()),
        reinterpret_cast<fftw_complex*>(x.data())
        );
}
