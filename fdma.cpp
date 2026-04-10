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
    // Генерация каждой поднесущей
    // -------------------------------------------------

    double sc_bb = 0;
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

        ScResult scRes;
        if(p.EnableSNRSym) scRes = sc.makeSc(symbolsPerCarrier[k].tr_sym_noisy, scp);
        else scRes = sc.makeSc(symbolsPerCarrier[k].tr_sym_clean, scp);

        maxLen = std::max(maxLen, scRes.tx.size());
        sc_bb = std::max(sc_bb, scRes.bandwidth);

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

    res.totalBandwidth =
        (p.FDMA_num_subcarriers - 1) * p.FDMA_step_carrier + sc_bb;

    res.perCarrierResults = carrierSignals;

    return res;
}

std::vector<std::vector<std::complex<double>>> FDMA::demodulate(
    const std::vector<std::complex<double>>& rxSignal, std::vector<Symbols>& symbols,
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

        auto scDemod = sc.demodulateSignal(rxSignal, symbols[k].tr_sym_clean, scp, scp);

        res.push_back(scDemod);
    }

    return res;
}

