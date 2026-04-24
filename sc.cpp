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
    res.bandwidth = Rs * (1.0 + p.SC_rolloff);

    // ----- Фильтр -----
    auto pulse = rrcFilter(span, sps, p.SC_rolloff);

    // ----- Polyphase filtering -----
    std::vector<std::complex<double>> baseband =
        polyphaseFilter(symbols, pulse, sps);

    // ----- Несущая -----
    res.tx.resize(baseband.size());
    res.t.resize(baseband.size());

    res.fc = p.fc;
    for(int n = 0; n < baseband.size(); ++n) {
        res.t[n] = double(n) / (p.fs * p.oversampling);

        double phase = 2.0 * M_PI * p.fc * res.t[n];

        res.tx[n] = baseband[n] *
                    std::exp(std::complex<double>(0.0, phase));
    }

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

std::vector<std::complex<double>> SC::demodulateSignal(
    const std::vector<std::complex<double>>& tx_signal, const std::vector<std::complex<double>> &sym_clean,
    const ScParams& p,
    ScParams& updatedParams)
{
    std::vector<std::complex<double>> baseband;

    double Rs = p.SC_symrate;
    int sps = std::round(p.fs * p.oversampling / Rs);

    // -------------------------------------------------
    // Downconversion
    // -------------------------------------------------
    baseband.resize(tx_signal.size());
    double invFs = 1.0 / (p.fs * p.oversampling);

    for (int n = 0; n < tx_signal.size(); ++n) {
        double t = n * invFs;

        baseband[n] = tx_signal[n] *
                      std::exp(std::complex<double>(0.0, -2.0 * M_PI * p.fc * t));
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

    if (!sym_clean.empty() &&
        sym_clean.size() == received.size())
    {
        std::complex<double> num(0.0, 0.0);
        std::complex<double> den(0.0, 0.0);

        for (int i = 0; i < received.size(); ++i) {
            num += std::conj(received[i]) * sym_clean[i];
            den += std::conj(received[i]) * received[i];
        }

        if (std::abs(den) > 1e-15) {
            std::complex<double> a = num / den;

            for (auto &s : received)
                s *= a;
        }
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
