#include "dpd.h"

DPD::DPD() {}

void DPD::train(const std::vector<std::complex<double> > &pa_input, const std::vector<std::complex<double> > &pa_output)
{
    int N = pa_output.size();
    int P = K * (M + 1);

    Eigen::MatrixXcd Phi(N, P);
    Eigen::VectorXcd x(N);

    for(int n=0; n<N; n++)
        x(n) = pa_input[n];

    for(int n=0; n<N; n++)
    {
        int col = 0;

        for(int m=0; m<=M; m++)
        {
            if(n-m < 0)
            {
                for(int k=0;k<K;k++)
                    Phi(n,col++) = 0;
                continue;
            }

            auto y = pa_output[n-m];
            double mag = std::abs(y);

            for(int k=1;k<=K;k++)
            {
                Phi(n,col++) = y * std::pow(mag, k-1);
            }
        }
    }

    Eigen::VectorXcd a =
        Phi.colPivHouseholderQr().solve(x);

    coeffs.resize(P);

    for(int i=0;i<P;i++)
        coeffs[i] = a(i);
}


std::vector<std::complex<double> > DPD::predistort(const std::vector<std::complex<double> > &input)
{
    int N = input.size();
    std::vector<std::complex<double>> output(N);

    for(int n=0;n<N;n++)
    {
        std::complex<double> sum = 0.0;

        int idx = 0;

        for(int m=0;m<=M;m++)
        {
            if(n-m < 0) continue;

            for(int k=1;k<=K;k++)
            {
                auto phi = basis(input[n-m],k);

                sum += coeffs[idx++] * phi;
            }
        }

        output[n] = sum;
    }

    return output;
}

void DPD::setOrder(int order)
{
    K = order;
}

void DPD::setMemory(int memory)
{
    M = memory;
}

std::complex<double> DPD::basis(const std::complex<double> &x, int k)
{
    double mag = std::abs(x);
    return x * std::pow(mag, k-1);
}
