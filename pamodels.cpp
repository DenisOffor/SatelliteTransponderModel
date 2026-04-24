#include "pamodels.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

PAModels::PAModels() {}

namespace
{
constexpr double EPS = 1e-15;

// В этой модели считаем, что dBm переводится в амплитуду как sqrt(mW):
// P_mW = 10^(P_dBm/10)
// A = sqrt(P_mW) = 10^(P_dBm/20)
double ampFromdBm(double p_dBm)
{
    return std::pow(10.0, p_dBm / 20.0);
}

double safe(double x)
{
    return std::max(x, EPS);
}

QString resolveStaticModel(const Source& source)
{
    if(source.PAModel == "Wiener")
        return source.W_StaticNonlinModel;

    if(source.PAModel == "Hammerstein")
        return source.H_StaticNonlinModel;

    if(source.PAModel == "Wiener-Hammerstein")
        return source.WH_StaticNonlinModel;

    return source.PAModel;
}

double calcAMAM(double r, const QString& model, const Source& source)
{
    if(model == "Saleh")
    {
        return source.SalehCoeffs[0] * r /
               (1.0 + source.SalehCoeffs[1] * r * r);
    }

    if(model == "Rapp")
    {
        return r /
               std::pow(
                   1.0 + std::pow(r / source.RappCoeffs[0],
                                  2.0 * source.RappCoeffs[1]),
                   1.0 / (2.0 * source.RappCoeffs[1])
                   );
    }

    if(model == "Ghorbani")
    {
        return source.GhorbaniCoeffs[0] * r /
               (
                   1.0
                   + source.GhorbaniCoeffs[1] * std::pow(r, 2.0)
                   + source.GhorbaniCoeffs[2] * std::pow(r, 4.0)
                   );
    }

    return r;
}

double calcAMPM(double r, const QString& model, const Source& source)
{
    if(model == "Saleh")
    {
        return source.SalehCoeffs[2] * r * r /
               (1.0 + source.SalehCoeffs[3] * r * r);
    }

    if(model == "Rapp")
    {
        return 0.0;
    }

    if(model == "Ghorbani")
    {
        return source.GhorbaniCoeffs[3] * r * r /
               (
                   1.0
                   + source.GhorbaniCoeffs[4] * std::pow(r, 2.0)
                   + source.GhorbaniCoeffs[5] * std::pow(r, 4.0)
                   );
    }

    return 0.0;
}

double findModelInputSat(const QString& model, const Source& source)
{
    if(model == "Saleh")
    {
        return 1.0 / std::sqrt(safe(source.SalehCoeffs[1]));
    }

    if(model == "Rapp")
    {
        // Для Rapp это параметр "knee/saturation input amplitude".
        return safe(source.RappCoeffs[0]);
    }

    if(model == "Ghorbani")
    {
        double rMax = 10.0;
        double step = 0.0005;

        double rSat = 0.0;
        double maxOut = -1.0;

        for(double r = 0.0; r <= rMax; r += step)
        {
            double y = calcAMAM(r, model, source);

            if(y > maxOut)
            {
                maxOut = y;
                rSat = r;
            }
        }

        return safe(rSat);
    }

    return 1.0;
}

void applyStaticModelWithPsat(std::vector<std::complex<double>>& sig,
                              const QString& staticModel,
                              Source& source)
{
    constexpr double EPS = 1e-15;

    auto dBmToAmp = [](double p_dBm) {
        // |x|^2 = mW
        // |x| = sqrt(mW)
        return std::pow(10.0, p_dBm / 20.0);
    };

    const double Pin_sat_amp  = std::max(dBmToAmp(source.Pin_sat_dBm), EPS);
    const double Pout_sat_amp = std::max(dBmToAmp(source.Pout_sat_dBm), EPS);

    /*
     * ВАЖНО:
     * r_ref_model — это точка модели, которую мы считаем "входным насыщением".
     * Именно ей соответствует Pin_sat_dBm.
     */
    double r_ref_model = 1.0;

    if(staticModel == "Saleh")
    {
        r_ref_model = 1.0 / std::sqrt(std::max(source.SalehCoeffs[1], EPS));
    }
    else if(staticModel == "Rapp")
    {
        r_ref_model = std::max(source.RappCoeffs[0], EPS);
    }
    else if(staticModel == "Ghorbani")
    {
        double rMax = 10.0;
        double step = 0.0005;

        double bestR = 0.0;
        double bestA = -1.0;

        for(double r = 0.0; r <= rMax; r += step)
        {
            double A =
                source.GhorbaniCoeffs[0] * r /
                (
                    1.0
                    + source.GhorbaniCoeffs[1] * std::pow(r, 2.0)
                    + source.GhorbaniCoeffs[2] * std::pow(r, 4.0)
                    );

            if(A > bestA)
            {
                bestA = A;
                bestR = r;
            }
        }

        r_ref_model = std::max(bestR, EPS);
    }

    auto calcAMAM = [&](double r) -> double
    {
        if(staticModel == "Saleh")
        {
            return source.SalehCoeffs[0] * r /
                   (1.0 + source.SalehCoeffs[1] * r * r);
        }
        else if(staticModel == "Rapp")
        {
            return r /
                   std::pow(
                       1.0 + std::pow(r / source.RappCoeffs[0],
                                      2.0 * source.RappCoeffs[1]),
                       1.0 / (2.0 * source.RappCoeffs[1])
                       );
        }
        else if(staticModel == "Ghorbani")
        {
            return source.GhorbaniCoeffs[0] * r /
                   (
                       1.0
                       + source.GhorbaniCoeffs[1] * std::pow(r, 2.0)
                       + source.GhorbaniCoeffs[2] * std::pow(r, 4.0)
                       );
        }

        return r;
    };

    auto calcAMPM = [&](double r) -> double
    {
        if(staticModel == "Saleh")
        {
            return source.SalehCoeffs[2] * r * r /
                   (1.0 + source.SalehCoeffs[3] * r * r);
        }
        else if(staticModel == "Rapp")
        {
            return 0.0;
        }
        else if(staticModel == "Ghorbani")
        {
            return source.GhorbaniCoeffs[3] * r * r /
                   (
                       1.0
                       + source.GhorbaniCoeffs[4] * std::pow(r, 2.0)
                       + source.GhorbaniCoeffs[5] * std::pow(r, 4.0)
                       );
        }

        return 0.0;
    };

    /*
     * Нормировочная выходная амплитуда модели.
     * Именно выход модели при r_ref_model будет соответствовать Pout_sat_dBm.
     */
    const double A_ref_model = std::max(calcAMAM(r_ref_model), EPS);

    for(size_t i = 0; i < sig.size(); ++i)
    {
        const double r_phys = std::abs(sig[i]);
        const double ph_in  = std::arg(sig[i]);

        /*
         * Если r_phys == Pin_sat_amp,
         * тогда r_model == r_ref_model.
         */
        const double r_model =
            r_phys / Pin_sat_amp * r_ref_model;

        const double A_model = calcAMAM(r_model);
        const double phi     = calcAMPM(r_model);

        /*
         * Если A_model == A_ref_model,
         * тогда A_out_phys == Pout_sat_amp.
         */
        const double A_out_phys =
            Pout_sat_amp * A_model / A_ref_model;

        sig[i] = std::polar(A_out_phys, ph_in + phi);
    }
}

} // namespace

double computeAvgPower(const std::vector<std::complex<double>>& x)
{
    if(x.empty())
        return 0.0;

    double sum = 0.0;

    for(const auto& s : x)
        sum += std::norm(s);

    return sum / static_cast<double>(x.size());
}

// ============================================================================
// Новые PA-модели с учетом Pin_sat_dBm / Pout_sat_dBm
// ============================================================================

void PAModels::SalehModel(std::vector<std::complex<double>>& sig,
                          Source& source)
{
    applyStaticModelWithPsat(sig, "Saleh", source);
}

void PAModels::RappModel(std::vector<std::complex<double>>& sig,
                         Source& source)
{
    applyStaticModelWithPsat(sig, "Rapp", source);
}

void PAModels::GhorbaniModel(std::vector<std::complex<double>>& sig,
                             Source& source)
{
    applyStaticModelWithPsat(sig, "Ghorbani", source);
}

void PAModels::WienerModel(std::vector<std::complex<double>>& sig,
                           QString Static_model,
                           std::vector<double>& FIR_Coeffs,
                           Source& source)
{
    ApplyFIRWithMemory(sig, FIR_Coeffs[0], 9);
    applyStaticModelWithPsat(sig, Static_model, source);
}

void PAModels::HammersteinModel(std::vector<std::complex<double>>& sig,
                                QString Static_model,
                                std::vector<double>& FIR_Coeffs,
                                Source& source)
{
    applyStaticModelWithPsat(sig, Static_model, source);
    ApplyFIRWithMemory(sig, FIR_Coeffs[0], 9);
}

void PAModels::WHModel(std::vector<std::complex<double>>& sig,
                       QString Static_model,
                       std::vector<double>& FIR_Coeffs,
                       Source& source)
{
    ApplyFIRWithMemory(sig, FIR_Coeffs[0], 9);
    applyStaticModelWithPsat(sig, Static_model, source);
    ApplyFIRWithMemory(sig, FIR_Coeffs[1], 9);
}

void PAModels::ApplyPA(std::vector<std::complex<double>>& sig,
                       Source& source)
{
    if(source.PAModel == "Saleh")
    {
        SalehModel(sig, source);
    }
    else if(source.PAModel == "Rapp")
    {
        RappModel(sig, source);
    }
    else if(source.PAModel == "Ghorbani")
    {
        GhorbaniModel(sig, source);
    }
    else if(source.PAModel == "Wiener")
    {
        WienerModel(sig,
                    source.W_StaticNonlinModel,
                    source.W_FIRCoeffs,
                    source);
    }
    else if(source.PAModel == "Hammerstein")
    {
        HammersteinModel(sig,
                         source.H_StaticNonlinModel,
                         source.H_FIRCoeffs,
                         source);
    }
    else if(source.PAModel == "Wiener-Hammerstein")
    {
        WHModel(sig,
                source.WH_StaticNonlinModel,
                source.WH_FIRCoeffs,
                source);
    }
}

// ============================================================================
// Вспомогательные функции
// ============================================================================

double PAModels::find_Asat_Ghorbani(const std::vector<double>& c,
                                    double gain_linear)
{
    Q_UNUSED(gain_linear);

    double Amax = 10.0;
    double step = 0.0005;

    double A_sat = 0.0;
    double max_out = -1.0;

    for(double Ain = 0.0; Ain < Amax; Ain += step)
    {
        double Aout =
            c[0] * Ain /
            (
                1.0
                + c[1] * Ain * Ain
                + c[2] * std::pow(Ain, 4.0)
                );

        if(Aout > max_out)
        {
            max_out = Aout;
            A_sat = Ain;
        }
    }

    return A_sat;
}

void PAModels::ApplyFIRWithMemory(std::vector<std::complex<double>>& signal,
                                  double alpha,
                                  int numTaps)
{
    const size_t N = signal.size();

    std::vector<std::complex<double>> h(numTaps);

    for(int m = 0; m < numTaps; ++m)
    {
        double mag = std::pow(alpha, m);
        double phase = alpha / 3.0 * m;
        h[m] = std::polar(mag, phase);
    }

    double sum = 0.0;

    for(int m = 0; m < numTaps; ++m)
        sum += std::abs(h[m]);

    if(sum > 0.0)
    {
        for(int m = 0; m < numTaps; ++m)
            h[m] /= sum;
    }

    std::vector<std::complex<double>> output(N, {0.0, 0.0});

    for(size_t n = 0; n < N; ++n)
    {
        for(int m = 0; m < numTaps; ++m)
        {
            if(n >= static_cast<size_t>(m))
                output[n] += h[m] * signal[n - m];
        }
    }

    signal = std::move(output);
}

void PAModels::ScaleToRMS_forPA(std::vector<std::complex<double>>& sig,
                                Source& source)
{
    if(sig.empty())
        return;

    // Теперь IBO считается относительно физического входного насыщения из GUI:
    //
    // Pin_work_dBm = Pin_sat_dBm - IBO_dB
    //
    // В амплитудах:
    // A_work_rms = A_sat / 10^(IBO_dB / 20)
    double Pin_sat_amp = ampFromdBm(source.Pin_sat_dBm);

    double target_rms =
        Pin_sat_amp / std::pow(10.0, source.IBO_dB / 20.0);

    scaleToRMS(sig, target_rms);
}

void PAModels::scaleToRMS(std::vector<std::complex<double>>& x,
                          double target_rms)
{
    if(x.empty())
        return;

    double sum = 0.0;

    for(const auto& s : x)
        sum += std::norm(s);

    double rms = std::sqrt(sum / static_cast<double>(x.size()));

    if(rms <= EPS)
        return;

    double k = target_rms / rms;

    for(auto& s : x)
        s *= k;
}
