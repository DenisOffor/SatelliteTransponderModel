#include "imuxomux.h"

#include <QFile>
#include <QTextStream>
#include <QStringList>
#include <QDebug>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace
{
constexpr double PI = 3.141592653589793238462643383279502884;
constexpr double EPS = 1e-15;

std::vector<std::complex<double>> ifftShiftEven(const std::vector<std::complex<double>>& x)
{
    const int N = static_cast<int>(x.size());
    std::vector<std::complex<double>> y(N);

    const int half = N / 2;
    for (int i = 0; i < half; ++i) {
        y[i] = x[half + i];
        y[half + i] = x[i];
    }
    return y;
}
}

ImuxOmux::ImuxOmux(const MuxFilterConfig& cfg)
{
    setConfig(cfg);
}

void ImuxOmux::setConfig(const MuxFilterConfig& cfg)
{
    m_cfg = cfg;

    m_imux = FilterState();
    m_omux = FilterState();

    m_imux.csvPath = cfg.imuxCsv;
    m_omux.csvPath = cfg.omuxCsv;
    m_imux.enabled = cfg.enableImux;
    m_omux.enabled = cfg.enableOmux;
}

bool ImuxOmux::recalc(double bandwidthHz, double sampleRateHz)
{
    if (bandwidthHz <= 0.0 || sampleRateHz <= 0.0) {
        qWarning() << "IMUX/OMUX: wrong bandwidth or sample rate" << bandwidthHz << sampleRateHz;
        return false;
    }

    if (m_cfg.synthesisNfft <= 0 || (m_cfg.synthesisNfft % 2) != 0) {
        qWarning() << "IMUX/OMUX: synthesisNfft must be positive and even" << m_cfg.synthesisNfft;
        return false;
    }

    if (m_cfg.taps <= 0 || m_cfg.taps >= m_cfg.synthesisNfft) {
        qWarning() << "IMUX/OMUX: taps must be in range 1..Nfft-1" << m_cfg.taps << m_cfg.synthesisNfft;
        return false;
    }

    bool ok = true;
    ok = recalcOne(m_imux, bandwidthHz, sampleRateHz) && ok;
    ok = recalcOne(m_omux, bandwidthHz, sampleRateHz) && ok;
    return ok;
}

void ImuxOmux::apply(std::vector<std::complex<double>>& signal,
                     MuxKind kind,
                     double bandwidthHz,
                     double sampleRateHz)
{
    if (signal.empty())
        return;

    recalc(bandwidthHz, sampleRateHz);

    const FilterState& s = state(kind);
    if (!s.enabled || !s.ready || s.h.empty())
        return;

    signal = applyFftFir(signal, s);

    if (m_cfg.compensateArtificialDelay)
        compensateDelay(signal);
}

const std::vector<std::complex<double>>& ImuxOmux::fir(MuxKind kind) const
{
    return state(kind).h;
}

bool ImuxOmux::isReady(MuxKind kind) const
{
    return state(kind).ready;
}

ImuxOmux::FilterState& ImuxOmux::state(MuxKind kind)
{
    return (kind == MuxKind::IMUX) ? m_imux : m_omux;
}

const ImuxOmux::FilterState& ImuxOmux::state(MuxKind kind) const
{
    return (kind == MuxKind::IMUX) ? m_imux : m_omux;
}

bool ImuxOmux::loadRaw(FilterState& s)
{
    if (!s.enabled)
        return false;

    if (s.rawLoaded)
        return true;

    QFile file(s.csvPath);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        qWarning() << "IMUX/OMUX: CSV not found, bypass for" << s.csvPath;
        s.ready = false;
        return false;
    }

    QTextStream in(&file);
    bool firstLine = true;

    struct Row { double f; double a; double gd; };
    std::vector<Row> rows;

    while (!in.atEnd()) {
        const QString line = in.readLine().trimmed();
        if (line.isEmpty())
            continue;

        if (firstLine) {
            firstLine = false;
            if (line.contains("freq", Qt::CaseInsensitive))
                continue;
        }

        const QStringList parts = line.split(',', Qt::SkipEmptyParts);
        if (parts.size() < 3)
            continue;

        bool okF = false;
        bool okA = false;
        bool okG = false;

        const double fMHz = parts[0].trimmed().toDouble(&okF);
        const double ampDb = parts[1].trimmed().toDouble(&okA);
        const double gdNs  = parts[2].trimmed().toDouble(&okG);

        if (okF && okA && okG)
            rows.push_back({fMHz * 1e6, ampDb, gdNs});
    }

    if (rows.size() < 2) {
        qWarning() << "IMUX/OMUX: CSV has not enough points" << s.csvPath;
        s.ready = false;
        return false;
    }

    std::sort(rows.begin(), rows.end(), [](const Row& l, const Row& r) {
        return l.f < r.f;
    });

    s.raw.fHz.clear();
    s.raw.ampDb.clear();
    s.raw.gdNs.clear();

    s.raw.fHz.reserve(rows.size());
    s.raw.ampDb.reserve(rows.size());
    s.raw.gdNs.reserve(rows.size());

    for (const Row& r : rows) {
        s.raw.fHz.push_back(r.f);
        s.raw.ampDb.push_back(r.a);
        s.raw.gdNs.push_back(r.gd);
    }

    s.rawLoaded = true;
    return true;
}

bool ImuxOmux::recalcOne(FilterState& s, double bandwidthHz, double sampleRateHz)
{
    if (!s.enabled) {
        s.ready = false;
        return true;
    }

    const bool same =
        s.ready &&
        almostEqual(s.lastBandwidthHz, bandwidthHz) &&
        almostEqual(s.lastSampleRateHz, sampleRateHz) &&
        s.lastSynthesisNfft == m_cfg.synthesisNfft &&
        s.lastTaps == m_cfg.taps &&
        almostEqual(s.lastReferenceBandwidthHz, m_cfg.referenceBandwidthHz);

    if (same)
        return true;

    if (!loadRaw(s))
        return false;

    s.h = synthesizeFir(s.raw, bandwidthHz, sampleRateHz);
    buildConvolutionSpectrum(s);

    s.lastBandwidthHz = bandwidthHz;
    s.lastSampleRateHz = sampleRateHz;
    s.lastSynthesisNfft = m_cfg.synthesisNfft;
    s.lastTaps = m_cfg.taps;
    s.lastReferenceBandwidthHz = m_cfg.referenceBandwidthHz;
    s.ready = !s.h.empty();

    return s.ready;
}

double ImuxOmux::interpLinearClipped(const std::vector<double>& x,
                                     const std::vector<double>& y,
                                     double xq)
{
    if (x.empty() || y.empty() || x.size() != y.size())
        return 0.0;

    if (xq <= x.front())
        return y.front();
    if (xq >= x.back())
        return y.back();

    auto it = std::lower_bound(x.begin(), x.end(), xq);
    const int i1 = static_cast<int>(std::distance(x.begin(), it));
    const int i0 = i1 - 1;

    const double dx = x[i1] - x[i0];
    if (std::abs(dx) < EPS)
        return y[i0];

    const double t = (xq - x[i0]) / dx;
    return y[i0] * (1.0 - t) + y[i1] * t;
}

void ImuxOmux::scaleResponse(const RawResponse& raw,
                             double referenceBandwidthHz,
                             double bandwidthHz,
                             const std::vector<double>& fGridHz,
                             std::vector<double>& ampDb,
                             std::vector<double>& gdNs)
{
    ampDb.resize(fGridHz.size());
    gdNs.resize(fGridHz.size());

    const double k = referenceBandwidthHz / bandwidthHz;
    const double fMin = raw.fHz.front();
    const double fMax = raw.fHz.back();
    const double ampFloor = *std::min_element(raw.ampDb.begin(), raw.ampDb.end());

    for (size_t i = 0; i < fGridHz.size(); ++i) {
        const double fQuery = fGridHz[i] * k;

        if (fQuery < fMin || fQuery > fMax) {
            ampDb[i] = ampFloor;
            gdNs[i] = 0.0;
            continue;
        }

        ampDb[i] = interpLinearClipped(raw.fHz, raw.ampDb, fQuery);
        gdNs[i] = k * interpLinearClipped(raw.fHz, raw.gdNs, fQuery);
    }
}

std::vector<std::complex<double>> ImuxOmux::synthesizeFir(const RawResponse& raw,
                                                          double bandwidthHz,
                                                          double sampleRateHz) const
{
    const int Nfft = m_cfg.synthesisNfft;
    const int Ntaps = m_cfg.taps;
    const int D = (Ntaps - 1) / 2;

    std::vector<double> f(Nfft);
    for (int i = 0; i < Nfft; ++i)
        f[i] = (i - Nfft / 2) * (sampleRateHz / static_cast<double>(Nfft));

    std::vector<double> ampDb;
    std::vector<double> gdNs;
    scaleResponse(raw, m_cfg.referenceBandwidthHz, bandwidthHz, f, ampDb, gdNs);

    std::vector<double> gdS(Nfft);
    for (int i = 0; i < Nfft; ++i)
        gdS[i] = gdNs[i] * 1e-9;

    // phi(f) = -2*pi * integral gd(f) df, аналог cumtrapz из Matlab.
    std::vector<double> phi(Nfft, 0.0);
    for (int i = 1; i < Nfft; ++i) {
        const double df = f[i] - f[i - 1];
        phi[i] = phi[i - 1] - 2.0 * PI * 0.5 * (gdS[i] + gdS[i - 1]) * df;
    }

    // phi(0)=0.
    int i0 = 0;
    double bestAbsF = std::abs(f[0]);
    for (int i = 1; i < Nfft; ++i) {
        const double af = std::abs(f[i]);
        if (af < bestAbsF) {
            bestAbsF = af;
            i0 = i;
        }
    }

    const double phi0 = phi[i0];
    for (double& p : phi)
        p -= phi0;

    // Делаем фазу нечётной: phi(-f) = -phi(f), как fliplr в Matlab-коде.
    std::vector<double> phiOdd(Nfft);
    for (int i = 0; i < Nfft; ++i)
        phiOdd[i] = 0.5 * (phi[i] - phi[Nfft - 1 - i]);
    phi.swap(phiOdd);

    // Добавляем искусственную постоянную задержку, чтобы FIR был причинным.
    for (int i = 0; i < Nfft; ++i)
        phi[i] -= 2.0 * PI * f[i] * (static_cast<double>(D) / sampleRateHz);

    // Комплексная ЧХ на fftshift-сетке.
    std::vector<std::complex<double>> Hshifted(Nfft);
    for (int i = 0; i < Nfft; ++i) {
        const double A = std::pow(10.0, ampDb[i] / 20.0);
        Hshifted[i] = std::polar(A, phi[i]);
    }

    // Matlab: hLong = ifft(ifftshift(H)).
    std::vector<std::complex<double>> hLong = ifftShiftEven(Hshifted);
    FFT fft(Nfft);
    fft.ifftInPlace(hLong);

    std::vector<std::complex<double>> h(Ntaps);
    for (int n = 0; n < Ntaps; ++n) {
        const double w = (Ntaps == 1)
            ? 1.0
            : 0.5 - 0.5 * std::cos(2.0 * PI * n / static_cast<double>(Ntaps - 1));
        h[n] = hLong[n] * w;
    }

    // Нормировка по DC: sum(h)=1.
    const std::complex<double> sumH = std::accumulate(
        h.begin(), h.end(), std::complex<double>(0.0, 0.0));

    if (std::abs(sumH) > EPS) {
        for (auto& v : h)
            v /= sumH;
    }

    return h;
}

void ImuxOmux::buildConvolutionSpectrum(FilterState& s) const
{
    const int Nfft = m_cfg.synthesisNfft;
    s.Hconv.assign(Nfft, std::complex<double>(0.0, 0.0));

    for (int i = 0; i < static_cast<int>(s.h.size()) && i < Nfft; ++i)
        s.Hconv[i] = s.h[i];

    FFT fft(Nfft);
    fft.fftInPlace(s.Hconv);
}

std::vector<std::complex<double>> ImuxOmux::applyFftFir(
    const std::vector<std::complex<double>>& x,
    const FilterState& s) const
{
    const int N = static_cast<int>(x.size());
    const int M = static_cast<int>(s.h.size());
    const int Nfft = m_cfg.synthesisNfft;

    if (N == 0 || M == 0)
        return x;

    if (M == 1) {
        std::vector<std::complex<double>> y = x;
        for (auto& v : y)
            v *= s.h[0];
        return y;
    }

    const int blockValid = Nfft - M + 1;
    if (blockValid <= 0)
        return x;

    std::vector<std::complex<double>> padded(M - 1 + N, std::complex<double>(0.0, 0.0));
    std::copy(x.begin(), x.end(), padded.begin() + (M - 1));

    std::vector<std::complex<double>> y(N, std::complex<double>(0.0, 0.0));
    std::vector<std::complex<double>> buf(Nfft);
    FFT fft(Nfft);

    int outPos = 0;
    while (outPos < N) {
        std::fill(buf.begin(), buf.end(), std::complex<double>(0.0, 0.0));

        for (int i = 0; i < Nfft; ++i) {
            const int idx = outPos + i;
            if (idx < static_cast<int>(padded.size()))
                buf[i] = padded[idx];
        }

        fft.fftInPlace(buf);
        for (int i = 0; i < Nfft; ++i)
            buf[i] *= s.Hconv[i];
        fft.ifftInPlace(buf);

        const int take = std::min(blockValid, N - outPos);
        for (int i = 0; i < take; ++i)
            y[outPos + i] = buf[M - 1 + i];

        outPos += take;
    }

    return y;
}

void ImuxOmux::compensateDelay(std::vector<std::complex<double>>& x) const
{
    const int D = (m_cfg.taps - 1) / 2;
    if (D <= 0 || x.empty())
        return;

    if (D >= static_cast<int>(x.size())) {
        std::fill(x.begin(), x.end(), std::complex<double>(0.0, 0.0));
        return;
    }

    for (int i = 0; i < static_cast<int>(x.size()) - D; ++i)
        x[i] = x[i + D];

    std::fill(x.end() - D, x.end(), std::complex<double>(0.0, 0.0));
}

bool ImuxOmux::almostEqual(double a, double b, double relTol)
{
    const double scale = std::max({1.0, std::abs(a), std::abs(b)});
    return std::abs(a - b) <= relTol * scale;
}
