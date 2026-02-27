#include "fdma.h"

FDMA::FDMA(SC& scRef)
    : sc(scRef) {}

FdmaResult FDMA::generate(
    const std::vector<Symbols> symbolsPerCarrier,
    const FdmaParams& p, const ScParams& sc_p)
{
    if(symbolsPerCarrier.size() != p.FDMA_num_subcarriers)
        throw std::runtime_error("Mismatch carriers count");

    FdmaResult res;

    std::vector<ScResult> carrierSignals;
    carrierSignals.reserve(p.FDMA_num_subcarriers);

    size_t maxLen = 0;

    // -------------------------------------------------
    // 1️⃣ Генерация каждой поднесущей
    // -------------------------------------------------

    for(int k = 0; k < p.FDMA_num_subcarriers; ++k)
    {
        ScParams scp;

        scp.SC_symrate = p.FDMA_symrate;
        scp.fs = p.fs;
        scp.oversampling = p.oversampling;
        scp.fc = p.FDMA_f_carrier + k * p.FDMA_step_carrier;
        scp.SNR_dB = 0;
        scp.SC_rolloff = sc_p.SC_rolloff;
        scp.SC_filter_length = sc_p.SC_filter_length;
        scp.SC_FilterType = sc_p.SC_FilterType;

        auto scRes = sc.makeSc(symbolsPerCarrier[k].tr_sym_noisy, scp);

        maxLen = std::max(maxLen, scRes.tx.size());

        carrierSignals.push_back(scRes);
    }

    // -------------------------------------------------
    // Выравнивание длины
    // -------------------------------------------------

    res.tx.resize(maxLen);
    std::fill(res.tx.begin(), res.tx.end(), std::complex<double>(0,0));

    for(const auto& sig : carrierSignals)
    {
        for(qsizetype n = 0; n < sig.tx.size(); ++n)
            res.tx[n] += sig.tx[n];
    }

    // -------------------------------------------------
    // Временная шкала
    // -------------------------------------------------

    res.t.resize(maxLen);

    double Fs_eff = static_cast<double>(p.fs) * p.oversampling;

    for(qsizetype n = 0; n < maxLen; ++n)
        res.t[n] = static_cast<double>(n) / Fs_eff;

    // -------------------------------------------------
    // Полоса
    // -------------------------------------------------

    if(p.SNRSig > 0)
        addAwgn(res, p.SNRSig);

    res.totalBandwidth =
        p.FDMA_num_subcarriers * p.FDMA_step_carrier;

    res.perCarrierResults = carrierSignals;

    return res;
}

std::vector<std::vector<std::complex<double>>> FDMA::demodulate(
    const std::vector<std::complex<double>>& rxSignal,
    const FdmaParams& p, const ScParams& sc_p)
{
    std::vector<std::vector<std::complex<double>>> res;
    res.reserve(p.FDMA_num_subcarriers);

    // -------------------------------------------------
    // Демодуляция каждой поднесущей
    // -------------------------------------------------

    for(int k = 0; k < p.FDMA_num_subcarriers; ++k)
    {
        ScParams scp;

        scp.SC_symrate = p.FDMA_symrate;
        scp.fs = p.fs;
        scp.oversampling = p.oversampling;
        scp.fc = p.FDMA_f_carrier + k * p.FDMA_step_carrier;
        scp.SNR_dB = p.SNRSig;
        scp.SC_rolloff = sc_p.SC_rolloff;
        scp.SC_filter_length = sc_p.SC_filter_length;
        scp.SC_FilterType = sc_p.SC_FilterType;

        auto scDemod = sc.demodulateSignal(rxSignal, scp, scp);

        res.push_back(scDemod);
    }

    return res;
}

void FDMA::addAwgn(FdmaResult &x, double SNR_dB)
{
    // Считаем мощность сигнала
    double power = 0;
    for(const auto& v : x.tx)
        power += std::norm(v);  // norm = real^2 + imag^2
    power /= x.tx.size();

    // Считаем дисперсию шума
    double snr = std::pow(10.0, SNR_dB / 10.0);
    double noiseVar = power / snr;
    double noiseStd = std::sqrt(noiseVar / 2);  // /2 потому что шум комплексный (I и Q компоненты)

    // Генератор шума
    static std::default_random_engine gen(std::random_device{}());  // static чтоб каждый раз не создавать
    std::normal_distribution<double> dist(0.0, noiseStd);

    // Добавляем шум
    x.currentNoise.resize(x.tx.size());
    for(int i = 0; i < x.tx.size(); ++i) {
        double noiseI = dist(gen);  // действительная часть
        double noiseQ = dist(gen);  // мнимая часть
        x.currentNoise[i] = std::complex<double>(noiseI, noiseQ);
        x.tx[i] += x.currentNoise[i];
    }
}

void FDMA::changeAwgn(FdmaResult &x, FdmaParams &p)
{
    for(int i = 0; i < x.tx.size(); ++i)
        x.tx[i] -= x.currentNoise[i];

    // Считаем мощность сигнала
    double power = 0;
    for(const auto& v : x.tx)
        power += std::norm(v);  // norm = real^2 + imag^2
    power /= x.tx.size();

    // Считаем дисперсию шума
    double snr = std::pow(10.0, p.SNRSig / 10.0);
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
