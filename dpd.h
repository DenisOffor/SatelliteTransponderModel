#ifndef DPD_H
#define DPD_H

#pragma once
#include <vector>
#include <complex>
#include <Eigen/Dense>

class DPD
{
public:

    DPD(int order, int memory)
        : K(order),
        M(memory),
        P(((order+1)/2)*(memory+1))
    {
        coeffs = Eigen::VectorXcd::Zero(P);
    }


    void train(const std::vector<std::complex<double>>& desired,
               const std::vector<std::complex<double>>& pa_output)
    {
        int N = desired.size();

        std::vector<std::complex<double>> pred = desired;

        double mu = 0.3;

        for(int iter=0; iter<5; iter++)
        {
            for(int n=0;n<N;n++)
            {
                std::complex<double> e = desired[n] - pa_output[n];
                pred[n] += mu * e;
            }
        }

        fit_polynomial(desired, pred);
    }

    void setMemory(int Memory) {
        this->M = Memory;
        P = K*(M+1);
        coeffs = Eigen::VectorXcd::Zero(P);
    }

    void setOrder(int order) {
        this->K = order;
        P = K*(M+1);
        coeffs = Eigen::VectorXcd::Zero(P);
    }

    std::vector<std::complex<double>> apply(
        const std::vector<std::complex<double>>& x)
    {
        int N = x.size();
        std::vector<std::complex<double>> y(N);

        for(int n=0;n<N;n++)
        {
            std::complex<double> acc=0;
            int c=0;

            for(int m=0;m<=M;m++)
            {
                if(n-m<0)
                {
                    c+=(K+1)/2;
                    continue;
                }

                auto v = x[n-m];
                double mag = std::norm(v);
                auto base = v;

                for(int k=1;k<=K;k+=2)
                {
                    acc += coeffs[c++] * base;
                    base *= mag;
                }
            }

            y[n]=acc;
        }

        return y;
    }


private:

    int K;
    int M;
    int P;

    Eigen::VectorXcd coeffs;


    void fit_polynomial(
        const std::vector<std::complex<double>>& x,
        const std::vector<std::complex<double>>& y)
    {
        int N = x.size();
        int rows = N-M;

        Eigen::MatrixXcd H(rows,P);
        Eigen::VectorXcd b(rows);

        for(int n=M;n<N;n++)
        {
            b[n-M] = y[n];

            int c=0;

            for(int m=0;m<=M;m++)
            {
                auto v = x[n-m];

                double mag = std::norm(v);
                auto base = v;

                for(int k=1;k<=K;k+=2)
                {
                    H(n-M,c++) = base;
                    base *= mag;
                }
            }
        }

        coeffs =
            (H.adjoint()*H).ldlt().solve(H.adjoint()*b);
    }
};

#endif // DPD_H
