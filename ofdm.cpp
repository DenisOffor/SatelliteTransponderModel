#include "ofdm.h"

OFDM::OFDM() {}

OfdmResult OFDM::makeOfdm(const QVector<std::complex<double> > &symbols, const OfdmParams &p)
{
    OfdmResult res;

    int Nfft_os = p.Nfft * p.oversampling;
    res.Tsym = double(Nfft_os) / p.fs;

    int Nactive = 0;
    auto X = ofdmSubcarrierMapping(symbols, p.Nfft,
                                   p.GB_DC, p.GB_Nyq, Nactive);
    res.Nactive = Nactive;

    // Oversampling (вставка нулей)
    if(p.oversampling > 1) {
        for(auto& row : X) {
            QVector<std::complex<double>> Xos(Nfft_os);
            int half = p.Nfft / 2;
            for(int i = 0; i < half; ++i)
                Xos[i] = row[i];
            for(int i = 0; i < half; ++i)
                Xos[Nfft_os - half + i] = row[half + i];
            row = Xos;
        }
    }

    res.BB = double(Nactive) / p.Nfft * p.fs;

    QVector<std::complex<double>> yBB;

    FFT myfft(X[0].size());
    for(const auto& row : X) {
        auto y = myfft.ifft(row);

        if(p.CP > 0) {
            QVector<std::complex<double>> ycp;
            for(int i = y.size() - p.CP; i < y.size(); ++i)
                ycp.append(y[i]);
            ycp += y;
            y = ycp;
        }
        yBB += y;
    }

    int N = yBB.size();
    res.t.resize(N);
    res.tx.resize(N);

    res.fc = p.fc;
    for(int n = 0; n < N; ++n) {
        res.t[n] = double(n) / (p.fs * p.oversampling);

        std::complex<double> s = yBB[n];
    //    if(p.fc > 0) {
    //        s *= std::exp(std::complex<double>(0, 2 * M_PI * p.fc * res.t[n]));
    //        res.tx[n] = std::real(s);
    //    } else {
            res.tx[n] = s;
    //    }
    }

    if(p.SNR_dB > 0)
        addAwgn(res, p.SNR_dB);

    return res;
}

QVector<std::complex<double>> OFDM::ofdm_demodulate(const QVector<std::complex<double>> &rx, const QVector<std::complex<double>> &sym_clean, const OfdmParams &p)
{
    int Nfft_os = p.Nfft * p.oversampling;
    int sym_len = Nfft_os + p.CP;

    int numSym = rx.size() / sym_len;

    QVector<std::complex<double>> rx_bb;
    rx_bb = rx;

    if (p.fc > 0) {
    //    rx_bb.resize(rx.size());
    //    double invFs = 1.0 / (p.fs * p.oversampling);

    //    for (int n = 0; n < rx.size(); ++n) {
    //        double t = n * invFs;
    //        std::complex<double> carrier =
    //            std::exp(std::complex<double>(0, -2 * M_PI * p.fc * t));

    //       rx_bb[n] = rx[n] * carrier;
    //    }
    }
    //else {
    //    rx_bb = rx;   // shallow copy (Qt implicit sharing)
    //}


    QVector<QVector<std::complex<double>>> symbols(numSym, QVector<std::complex<double>>(Nfft_os));

    // Разбивка + удаление CP
    for (int s = 0; s < numSym; ++s)
    {
        int offset = s * sym_len + p.CP;

        for (int n = 0; n < Nfft_os; ++n)
            symbols[s][n] = rx[offset + n];
    }

    FFT myfft(symbols[0].size());
    // FFT по каждому символу
    for (int s = 0; s < numSym; ++s) {
        symbols[s] = myfft.fft(symbols[s]);
    }

    // Убираем oversampling (если есть)
    if (p.oversampling > 1)
    {
        int half = p.Nfft / 2;

        for (int s = 0; s < numSym; ++s)
        {
            QVector<std::complex<double>> tmp(p.Nfft);

            // нижняя половина
            for (int k = 0; k < half; ++k)
                tmp[k] = symbols[s][k];

            // верхняя половина
            for (int k = 0; k < half; ++k)
                tmp[half + k] =
                    symbols[s][Nfft_os - half + k];

            symbols[s] = tmp;
        }
    }
    QVector<std::complex<double>> symbols_rx =
        ofdm_subcarrier_demapping(symbols, p.Nfft, p.GB_DC, p.GB_Nyq, sym_clean.size());

    if (!sym_clean.isEmpty() &&
        sym_clean.size() == symbols_rx.size())
    {
        std::complex<double> num(0.0, 0.0);
        std::complex<double> den(0.0, 0.0);

        for (int i = 0; i < symbols_rx.size(); ++i) {
            num += std::conj(symbols_rx[i]) * sym_clean[i];
            den += std::conj(symbols_rx[i]) * symbols_rx[i];
        }

        if (std::abs(den) > 1e-15) {
            std::complex<double> a = num / den;

            for (auto &s : symbols_rx)
                s *= a;
        }
    }

    return symbols_rx;
}

QVector<QVector<std::complex<double>>> OFDM::ofdmSubcarrierMapping(const QVector<std::complex<double> > &dataSymbols, const int Nfft,
                                                                   const int GB_DC, const int GB_Nyq, int &Nactive)
{
    int half = Nfft / 2;

    int pos_start = GB_DC;
    int pos_end   = half - GB_Nyq - 1;

    int neg_start = half + GB_Nyq;
    int neg_end   = Nfft - GB_DC - 1;

    QVector<int> activeIdx;
    for(int i = pos_start; i <= pos_end; ++i) activeIdx.append(i);
    for(int i = neg_start; i <= neg_end; ++i) activeIdx.append(i);

    Nactive = activeIdx.size();

    int numRows = std::ceil(double(dataSymbols.size()) / Nactive);

    QVector<std::complex<double>> padded = dataSymbols;
    padded.resize(numRows * Nactive); // zero-padding

    QVector<QVector<std::complex<double>>> X(numRows,
                                             QVector<std::complex<double>>(Nfft, {0.0, 0.0}));

    int k = 0;
    for(int r = 0; r < numRows; ++r)
        for(int c : activeIdx)
            X[r][c] = padded[k++];

    return X;
}

QVector<std::complex<double> > OFDM::ofdm_subcarrier_demapping(const QVector<QVector<std::complex<double> > > &X, const int Nfft, const int GB_DC,
                                                              const int GB_Nyq, const int NumSym)
{
    int half = Nfft / 2;

    int pos_start = GB_DC;
    int pos_end   = half - GB_Nyq - 1;

    int neg_start = half + GB_Nyq;
    int neg_end   = Nfft - GB_DC - 1;

    QVector<std::complex<double>> data;

    for (const auto& row : X)
    {
        // положительные частоты
        for (int k = pos_start; k <= pos_end; ++k) {
            if (data.size() >= NumSym)  // проверяем ДО добавления
                return data;
            data.push_back(row[k]);
        }

        // отрицательные частоты
        for (int k = neg_start; k <= neg_end; ++k) {
            if (data.size() >= NumSym)  // проверяем ДО добавления
                return data;
            data.push_back(row[k]);
        }
    }

    return data;
}

void OFDM::addAwgn(OfdmResult &x, double SNR_dB)
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

void OFDM::changeFc(OfdmResult &x, OfdmParams& p)
{
    double invFs = 1.0 / (p.fs * p.oversampling);

    // Вектор для хранения комплексной огибающей
    std::vector<std::complex<double>> baseband(x.tx.size());

    for (int n = 0; n < x.tx.size(); ++n) {
        double t = n * invFs;

        if(x.fc > 0) {
            // Шаг 1: Демодуляция в baseband (убираем старую несущую)
            std::complex<double> mixer_old =
                std::exp(std::complex<double>(0, -2 * M_PI * x.fc * t));
            baseband[n] = x.tx[n] * mixer_old;
        } else {
            // Если старой несущей нет, сигнал уже в baseband
            baseband[n] = x.tx[n];
        }
    }

    for (int n = 0; n < x.tx.size(); ++n) {
        double t = n * invFs;

        if(p.fc > 0) {
            // Шаг 2: Модуляция на новую несущую
            std::complex<double> mixer_new =
                std::exp(std::complex<double>(0, 2 * M_PI * p.fc * t));
            x.tx[n] = std::real(baseband[n] * mixer_new);
        } else {
            // Если новой несущей нет, оставляем в baseband
            x.tx[n] = baseband[n];
        }
    }

    x.fc = p.fc;
}

void OFDM::changeAwgn(OfdmResult &x, OfdmParams &p)
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
