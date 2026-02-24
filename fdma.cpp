#include "fdma.h"

FDMA::FDMA(SC& scRef)
    : sc(scRef) {}

FdmaResult FDMA::generate(
    const QVector<Symbols> symbolsPerCarrier,
    const FdmaParams& p)
{
    if(symbolsPerCarrier.size() != p.FDMA_num_subcarriers)
        throw std::runtime_error("Mismatch carriers count");

    FdmaResult res;

    QVector<ScResult> carrierSignals;
    carrierSignals.reserve(p.FDMA_num_subcarriers);

    qsizetype maxLen = 0;

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
        scp.SNR_dB = p.SNRSig;
        scp.SC_rolloff = p.rolloff;
        scp.SC_filter_length = p.filterSpan;
        scp.SC_FilterType = p.filterType;

        auto scRes = sc.makeSc(symbolsPerCarrier[k].tr_sym_noisy, scp);

        // ✅ FIX: теперь оба аргумента одного типа
        maxLen = std::max(maxLen, scRes.tx.size());

        carrierSignals.append(scRes);
    }

    // -------------------------------------------------
    // 2️⃣ Выравнивание длины
    // -------------------------------------------------

    res.tx.resize(maxLen);
    std::fill(res.tx.begin(), res.tx.end(), std::complex<double>(0,0));

    for(const auto& sig : carrierSignals)
    {
        for(qsizetype n = 0; n < sig.tx.size(); ++n)
            res.tx[n] += sig.tx[n];
    }

    // -------------------------------------------------
    // 3️⃣ Временная шкала
    // -------------------------------------------------

    res.t.resize(maxLen);

    double Fs_eff = static_cast<double>(p.fs) * p.oversampling;

    for(qsizetype n = 0; n < maxLen; ++n)
        res.t[n] = static_cast<double>(n) / Fs_eff;

    // -------------------------------------------------
    // 4️⃣ Полоса
    // -------------------------------------------------

    res.totalBandwidth =
        p.FDMA_num_subcarriers * p.FDMA_step_carrier;

    res.perCarrierResults = carrierSignals;

    return res;
}
