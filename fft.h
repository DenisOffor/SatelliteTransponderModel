#ifndef FFT_H
#define FFT_H
#include <complex>
#include <cmath>
#include <fftw3.h>
#include "HelpfullStructs.h"

class FFT
{
public:
    FFT(int n);
    ~FFT();

    //std::vector<std::complex<double>> fft(const std::vector<std::complex<double>>& x);
    //std::vector<std::complex<double>> ifft(const std::vector<std::complex<double>>& X);
    void fftInPlace(std::vector<std::complex<double>> x);
    void ifftInPlace(std::vector<std::complex<double>> x);

private:
    int m_n;
    std::vector<std::complex<double>> m_buffer;
    fftw_plan m_planForward;
    fftw_plan m_planBackward;
    bool m_plansCreated;
};


#endif // FFT_H
