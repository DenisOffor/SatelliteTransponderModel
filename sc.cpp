#include "sc.h"

SC::SC() {}

ScResult SC::makeSc(const std::vector<std::complex<double>>& inputSymbols,
                    const ScParams& p)
{
    ScResult res;

    double Rs = p.SC_symrate;
    int sps = std::round(p.fs * p.oversampling / Rs);

    if (sps < 2)
        throw std::runtime_error("Fs too low for symbol rate");

    // ----- Паддинг символов -----
    std::vector<std::complex<double>> symbols;

    int span = p.SC_filter_length;

    // padding в начале
    for(int i = 0; i < span; ++i)
        symbols.push_back(std::complex<double>(0,0));

    symbols.insert(symbols.end(), inputSymbols.begin(), inputSymbols.end());

    for(int i = 0; i < span; ++i)
        symbols.push_back(std::complex<double>(0,0));

    // ----- Полоса -----
    res.bandwidth = (p.SC_FilterType.toLower() == "rrc") ?
                        Rs * (1.0 + p.SC_rolloff) :
                        Rs;

    // ----- Фильтр -----
    auto pulse = rrcFilter(span, sps, p.SC_rolloff);

    // ----- Polyphase filtering -----
    std::vector<std::complex<double>> baseband =
        polyphaseFilter(symbols, pulse, sps);

    // ----- Несущая -----
    res.tx.resize(baseband.size());
    res.t.resize(baseband.size());

    double Fs_eff = p.fs * p.oversampling;

    res.fc = p.fc;

    // Fast carrier generation (recursive phase rotation)
    double w = 2 * M_PI * p.fc / Fs_eff;

    std::complex<double> carrier_step =
        std::exp(std::complex<double>(0, w));

    std::complex<double> carrier = 1.0;

    for(int n=0;n<baseband.size();n++)
    {
        double t = double(n) / Fs_eff;
        res.t[n] = t;

        if(p.fc > 0)
        {
            std::complex<double> mod =
                baseband[n] * carrier;

            res.tx[n] = std::real(mod);

            carrier *= carrier_step;
        }
        else
        {
            res.tx[n] = baseband[n];
        }
    }

    if(p.SNR_dB > 0)
        addAwgn(res, p.SNR_dB);

    return res;
}

std::vector<std::complex<double>> SC::polyphaseFilter(
    const std::vector<std::complex<double>>& symbols,
    const std::vector<double>& h,
    int sps)
{
    int L = h.size();
    int phases = sps;

    std::vector<std::vector<double>> poly(phases);

    // Polyphase decomposition
    for(int p=0;p<phases;p++)
        for(int k=p;k<L;k+=phases)
            poly[p].push_back(h[k]);

    int delay = L/2;

    std::vector<std::complex<double>> out(symbols.size()*sps);
    std::fill(out.begin(), out.end(), std::complex<double>(0,0));

    int N = symbols.size();

    for(int n=0;n<N;n++)
    {
        for(int p=0;p<phases;p++)
        {
            std::complex<double> acc = 0;

            auto &phaseFilter = poly[p];

            for(int k=0;k<phaseFilter.size();k++)
            {
                int sym_idx = n - k;

                if(sym_idx >= 0)
                    acc += symbols[sym_idx] * phaseFilter[k];
            }

            int idx = n*sps + p;

            if(idx < out.size())
                out[idx] = acc;
        }
    }

    return out;
}

std::vector<std::complex<double>> SC::upsample(
    const std::vector<std::complex<double>>& in, int sps)
{
    std::vector<std::complex<double>> out(in.size() * sps);

    for(int i = 0; i < in.size(); ++i)
        out[i*sps] = in[i];

    return out;
}

std::vector<double> SC::rrcFilter(int span, int sps, double beta)
{
    int N = span * sps;
    std::vector<double> h(N+1);

    double T = 1.0;
    int mid = N/2;

    for(int n = 0; n <= N; ++n)
    {
        double t = (n - mid) / double(sps);

        if (std::abs(t) < 1e-8)
        {
            h[n] = (1.0 + beta*(4/M_PI - 1));
        }
        else if (std::abs(std::abs(t) - T/(4*beta)) < 1e-8)
        {
            double val = (beta/std::sqrt(2)) *
                         ((1+2/M_PI)*std::sin(M_PI/(4*beta)) +
                          (1-2/M_PI)*std::cos(M_PI/(4*beta)));
            h[n] = val;
        }
        else
        {
            double num = std::sin(M_PI*t*(1-beta)/T)
            + 4*beta*t/T*std::cos(M_PI*t*(1+beta)/T);
            double den = M_PI*t*(1 - std::pow(4*beta*t/T,2))/T;
            h[n] = num / den;
        }
    }

    // нормировка энергии = 1
    double energy = 0;
    for(double v : h) energy += v*v;
    energy = std::sqrt(energy);
    for(auto &v : h) v /= energy;

    return h;
}

std::vector<std::complex<double>> SC::filterSignal(
    const std::vector<std::complex<double>>& x,
    const std::vector<double>& h)
{
    std::vector<std::complex<double>> y(x.size());

    for(int n = 0; n < x.size(); ++n)
    {
        std::complex<double> acc = 0;
        for(int k = 0; k < h.size(); ++k)
        {
            if(n-k >= 0)
                acc += x[n-k] * h[k];
        }
        y[n] = acc;
    }

    return y;
}

void SC::addAwgn(ScResult &x, double SNR_dB)
{
    double power = 0;
    for(const auto& v : x.tx)
        power += std::norm(v);
    power /= x.tx.size();

    double snr = std::pow(10.0, SNR_dB/10.0);
    double noiseVar = power / snr;
    double noiseStd = std::sqrt(noiseVar/2);

    static std::default_random_engine gen(std::random_device{}());
    std::normal_distribution<double> dist(0.0, noiseStd);

    x.currentNoise.resize(x.tx.size());

    for(int i=0;i<x.tx.size();++i)
    {
        double ni = dist(gen);
        double nq = dist(gen);

        x.currentNoise[i] = {ni, nq};
        x.tx[i] += x.currentNoise[i];
    }
}

std::vector<std::complex<double>> SC::demodulateSignal(
    const std::vector<std::complex<double>>& tx_signal,
    const ScParams& p,
    ScParams& updatedParams)
{
    std::vector<std::complex<double>> baseband;

    double Rs = p.SC_symrate;
    int sps = std::round(p.fs * p.oversampling / Rs);
    double Fs_eff = p.fs * p.oversampling;

    // -------------------------------------------------
    // 1️⃣ Downconversion
    // -------------------------------------------------

    if (p.fc > 0)
    {
        baseband.resize(tx_signal.size());

        double w = 2 * M_PI * p.fc / Fs_eff;
        std::complex<double> step =
            std::exp(std::complex<double>(0, -w)); // минус для downconversion

        std::complex<double> carrier = 1.0;

        for(int n = 0; n < tx_signal.size(); ++n)
        {
            baseband[n] = tx_signal[n] * carrier;
            carrier *= step;
        }
    }
    else
    {
        baseband = tx_signal;
    }

    // -------------------------------------------------
    // Matched filter (тот же RRC)
    // -------------------------------------------------

    int span = p.SC_filter_length;

    std::vector<double> rx_filter;

    if (p.SC_FilterType.toLower() == "rrc")
        rx_filter = rrcFilter(span, sps, p.SC_rolloff);
    else
    {
        rx_filter = std::vector<double>(sps, 1.0);
        double norm = 0;
        for(double v : rx_filter) norm += v*v;
        norm = std::sqrt(norm);
        for(auto &v : rx_filter) v /= norm;
    }

    auto received = matchedFilterDecimate(baseband, rx_filter, sps);

    // удалить суммарную задержку (span символов)
    int total_delay = 2*span;

    if(received.size() > total_delay)
    {
        received.erase(received.begin(),
                       received.begin() + total_delay);
    }


    // -------------------------------------------------
    // Нормализация мощности
    // -------------------------------------------------

    double power = 0;

    for(const auto& v : received)
        power += std::norm(v);

    if(received.size() > 0)
        power /= received.size();

    if(power > 0)
    {
        double scale = std::sqrt(power);

        for(auto& v : received)
            v /= scale;
    }

    return received;
}

std::vector<std::complex<double>> SC::matchedFilterDecimate(
    const std::vector<std::complex<double>>& x,
    const std::vector<double>& h,
    int sps)
{
    int L = h.size();
    int N = x.size();

    int outSize = N / sps;
    std::vector<std::complex<double>> out;
    out.reserve(outSize);

    int delay = L / 2;

    for(int n = 0; n < N; n += sps)
    {
        std::complex<double> acc = 0;

        for(int k = 0; k < L; ++k)
        {
            int idx = n - k;
            if(idx >= 0)
                acc += x[idx] * h[k];
        }

        out.push_back(acc);
    }

    return out;
}

void SC::changeAwgn(ScResult &x, ScParams &p)
{
    for(int i = 0; i < x.tx.size(); ++i)
        x.tx[i] -= x.currentNoise[i];

    // Считаем мощность сигнала
    double power = 0;
    for(const auto& v : x.tx)
        power += std::norm(v);  // norm = real^2 + imag^2
    power /= x.tx.size();

    // Считаем дисперсию шума
    double snr = std::pow(10.0, p.SNR_dB / 10.0);
    double noiseVar = power / snr;
    double noiseStd = std::sqrt(noiseVar / 2);  // /2 потому что шум комплексный (I и Q компоненты)

    // Генератор шума
    static std::default_random_engine gen(std::random_device{}());  // static чтоб каждый раз не создавать
    std::normal_distribution<double> dist(0.0, noiseStd);

    // Добавляем шум

    for(int i = 0; i < x.tx.size(); ++i) {
        double noiseI = dist(gen);  // действительная часть
        double noiseQ = dist(gen);  // мнимая часть
        x.currentNoise[i] = std::complex<double>(noiseI, noiseQ);
        x.tx[i] += x.currentNoise[i];
    }
}
