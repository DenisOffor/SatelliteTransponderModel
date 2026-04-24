#ifndef IMUXOMUX_H
#define IMUXOMUX_H

#include <complex>
#include <vector>
#include <QString>
#include <fft.h>
#include <QCoreApplication>
#include <QDir>

// Какой именно фильтр применяем.
enum class MuxKind
{
    IMUX,
    OMUX
};

struct MuxFilterConfig
{
    // Пути к CSV с колонками: freq_MHz, amp_dB, group_delay_ns.
    // Если файл IMUX пока отсутствует, IMUX будет работать как bypass.
    QString filtersDir =
        QCoreApplication::applicationDirPath() + "/filters";

    QString imuxCsv = filtersDir + "/imux_digitized_approx.csv";
    QString omuxCsv = filtersDir + "/omux_digitized_approx.csv";

    // Исходная характеристика оцифрована для 36 МГц.
    double referenceBandwidthHz = 36.0e6;

    // Параметры синтеза FIR, как в Matlab-файле.
    int synthesisNfft = 8192;
    int taps = 51;

    // Для твоей модели лучше убирать только искусственную задержку D=(Ntaps-1)/2,
    // иначе демодуляторы SC/OFDM/FDMA получат сдвиг по времени.
    bool compensateArtificialDelay = true;

    bool enableImux = true;
    bool enableOmux = true;
};

class ImuxOmux
{
public:
    explicit ImuxOmux(const MuxFilterConfig& cfg = MuxFilterConfig());

    void setConfig(const MuxFilterConfig& cfg);

    // Внешний вызов №1: пересчитать оба FIR под новую полосу.
    // Если bandwidthHz/Fs/параметры не изменились, ничего не делает.
    bool recalc(double bandwidthHz, double sampleRateHz);

    // Внешний вызов №2: применить IMUX или OMUX.
    // Внутри есть защита от лишнего пересчета при той же полосе.
    void apply(std::vector<std::complex<double>>& signal,
               MuxKind kind,
               double bandwidthHz,
               double sampleRateHz);

    const std::vector<std::complex<double>>& fir(MuxKind kind) const;
    bool isReady(MuxKind kind) const;

private:
    struct RawResponse
    {
        std::vector<double> fHz;
        std::vector<double> ampDb;
        std::vector<double> gdNs;
    };

    struct FilterState
    {
        QString csvPath;
        bool enabled = true;
        bool rawLoaded = false;
        bool ready = false;

        RawResponse raw;
        std::vector<std::complex<double>> h;
        std::vector<std::complex<double>> Hconv;

        double lastBandwidthHz = -1.0;
        double lastSampleRateHz = -1.0;
        int lastSynthesisNfft = -1;
        int lastTaps = -1;
        double lastReferenceBandwidthHz = -1.0;
    };

    MuxFilterConfig m_cfg;
    FilterState m_imux;
    FilterState m_omux;

    FilterState& state(MuxKind kind);
    const FilterState& state(MuxKind kind) const;

    bool loadRaw(FilterState& s);
    bool recalcOne(FilterState& s, double bandwidthHz, double sampleRateHz);

    static double interpLinearClipped(const std::vector<double>& x,
                                      const std::vector<double>& y,
                                      double xq);

    static void scaleResponse(const RawResponse& raw,
                              double referenceBandwidthHz,
                              double bandwidthHz,
                              const std::vector<double>& fGridHz,
                              std::vector<double>& ampDb,
                              std::vector<double>& gdNs);

    std::vector<std::complex<double>> synthesizeFir(const RawResponse& raw,
                                                    double bandwidthHz,
                                                    double sampleRateHz) const;

    void buildConvolutionSpectrum(FilterState& s) const;

    std::vector<std::complex<double>> applyFftFir(
        const std::vector<std::complex<double>>& x,
        const FilterState& s) const;

    void compensateDelay(std::vector<std::complex<double>>& x) const;

    static bool almostEqual(double a, double b, double relTol = 1e-12);
};

#endif // IMUXOMUX_H
