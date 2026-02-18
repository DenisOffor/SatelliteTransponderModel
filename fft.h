#ifndef FFT_H
#define FFT_H
#include <QVector>
#include <complex>
#include <cmath>
#include <fftw3.h>

class FFT
{
public:
    FFT(int n);
    ~FFT();

    QVector<std::complex<double>> fft(const QVector<std::complex<double>>& x);
    QVector<std::complex<double>> ifft(const QVector<std::complex<double>>& X);

private:
    int m_n;
    fftw_plan m_planForward;
    fftw_plan m_planBackward;
    fftw_complex *m_in, *m_out;
    bool m_plansCreated;
};


#endif // FFT_H
