#include "signalprocessing.h"

SignalProcessing::SignalProcessing() : myFdma(mySC), mydpd(), MyMux() {
    PACurve = new PaCurve(200);
    // Input power in dB and linear
    double dBm_start = -25;
    double dBm_end = 5;
    double dBm_step = (dBm_end - dBm_start) / (PACurve->point_num - 1);
    for(int i = 0; i < PACurve->point_num; i ++) {
        PACurve->P_in_abs.dB[i] = dBm_start + i * dBm_step;
        PACurve->P_in_abs.linear[i] = qPow(10.0, PACurve->P_in_abs.dB[i] / 10.0);
        PACurve->r_in[i] = qSqrt(PACurve->P_in_abs.linear[i]);
    }
    freq.resize(6);
    PSDs.resize(6);
    InitializeConstellations();
}


SignalProcessing::~SignalProcessing() {
    delete PACurve;
}

void SignalProcessing::MainLogicWork(NeedToRecalc CurrentRecalcNeeds)
{
    QElapsedTimer timer;
    timer.start();
    if(CurrentRecalcNeeds.DPDRecalc)
        RecalcDPD(CurrentRecalcNeeds);
    qDebug() << "DPD train time:" << timer.elapsed() << "ms";

    timer.restart();
    if(CurrentRecalcNeeds.RecalcSymbols || CurrentRecalcNeeds.RecalcNoiseSym || CurrentRecalcNeeds.FullRecalc)
        GeneratePacksOfSymbols(MySymbols, MySource, CurrentRecalcNeeds);
    qDebug() << "Sym time:" << timer.elapsed() << "ms";

    timer.restart();
    if(CurrentRecalcNeeds.RecalcSig || CurrentRecalcNeeds.FullRecalc)
        TransmitSignalProcessing(MySource, MySymbols, CurrentRecalcNeeds, CurrentRes);
    qDebug() << "Sig time:" << timer.elapsed() << "ms";

    timer.restart();
    if(CurrentRecalcNeeds.PARecalc || CurrentRecalcNeeds.FullRecalc || CurrentRecalcNeeds.RecalcNoiseSig)
            PAProcessing(MySource, CurrentRecalcNeeds, CurrentRes);
    qDebug() << "PA time:" << timer.elapsed() << "ms";

    timer.restart();
    if(CurrentRecalcNeeds.RecRecalc)
        ReceiveSignalProcessing(MySource, MySymbols, CurrentRecalcNeeds, CurrentRes);
    qDebug() << "Rec time:" << timer.elapsed() << "ms";

    timer.restart();
    if(CurrentRecalcNeeds.MetricsRecalc) {
        MyMetricsEval.computePSD(CurrentRes.tx_sig, MySource.fs, MySource.oversampling, freq[0], PSDs[0], OFDM_tx_sig_buff);
        MyMetricsEval.computePSD(CurrentRes.tx_plus_dpd_sig, MySource.fs, MySource.oversampling, freq[0], PSDs[1], OFDM_tx_plus_DPD_sig_buff);
        MyMetricsEval.computePSD(CurrentRes.pa_sig_noisy, MySource.fs, MySource.oversampling, freq[0], PSDs[2], OFDM_pa_sig_buff);
        MyMetricsEval.computePSD(CurrentRes.pa_plus_dpd_sig_noisy, MySource.fs, MySource.oversampling, freq[0], PSDs[3], OFDM_pa_plus_DPD_sig_buff);

        static int num_iter = 1;
        if(CurrentRecalcNeeds.CycleMode == true) {
            double BER_noDPD, BER_withDPD, EVM_noDPD, EVM_withDPD, Gain_noDPD, Gain_withDPD, NMSE_noDPD, NMSE_withDPD;
            double P_formed_noDPD, P_formed_withDPD, P_emitted_noDPD, P_emitted_withDPD, PARP_noDPD, PARP_withDPD;
            QPair<double, double> ACLR_noDPD, ACLR_withDPD;
            double CurrentOBO_noDPD, CurrentOBO_withDPD;
            std::tie(BER_noDPD, BER_withDPD) = MyMetricsEval.Calc_BER(MySymbols);
            std::tie(EVM_noDPD, EVM_withDPD) = MyMetricsEval.Calc_EVM(MySymbols);
            P_formed_noDPD = MyMetricsEval.compute_av_P_dBm(CurrentRes.tx_sig);
            P_formed_withDPD = MyMetricsEval.compute_av_P_dBm(CurrentRes.tx_plus_dpd_sig);
            P_emitted_noDPD = MyMetricsEval.compute_av_P_dBm(CurrentRes.pa_sig);
            P_emitted_withDPD = MyMetricsEval.compute_av_P_dBm(CurrentRes.pa_plus_dpd_sig);
            Gain_noDPD = CurrentRes.P_emitted_noDPD - CurrentRes.P_formed_noDPD;
            Gain_withDPD = CurrentRes.P_emitted_withDPD - CurrentRes.P_formed_withDPD;
            PARP_noDPD = MyMetricsEval.computePAPR_dB(CurrentRes.tx_sig);
            PARP_withDPD = MyMetricsEval.computePAPR_dB(CurrentRes.tx_plus_dpd_sig);
            NMSE_noDPD = MyMetricsEval.computeNMSE_dB(CurrentRes.tx_sig, CurrentRes.pa_sig);
            NMSE_withDPD = MyMetricsEval.computeNMSE_dB(CurrentRes.tx_sig, CurrentRes.pa_plus_dpd_sig);
            CurrentOBO_noDPD = MyMetricsEval.computeOBO_dB(CurrentRes.pa_sig, MySource.Pout_sat_dBm);
            CurrentOBO_withDPD = MyMetricsEval.computeOBO_dB(CurrentRes.pa_plus_dpd_sig_noisy, MySource.Pout_sat_dBm);

            ACLR_noDPD = MyMetricsEval.computeACPR(freq[0], PSDs[2], CurrentRes.BB, CurrentRes.BB * MySource.bb_delta, MySource.SigType);
            ACLR_withDPD = MyMetricsEval.computeACPR(freq[0], PSDs[3], CurrentRes.BB, CurrentRes.BB * MySource.bb_delta, MySource.SigType);

            CurrentRes.BER_noDPD = (CurrentRes.BER_noDPD * MySource.NumSym + BER_noDPD * MySource.NumSym) / ((++num_iter) * MySource.NumSym);
            CurrentRes.BER_withDPD = (CurrentRes.BER_withDPD * MySource.NumSym + BER_withDPD * MySource.NumSym) / (num_iter * MySource.NumSym);;

            CurrentRes.EVM_noDPD = (CurrentRes.EVM_noDPD * (num_iter - 1) + EVM_noDPD) / num_iter;
            CurrentRes.EVM_withDPD = (CurrentRes.EVM_withDPD * (num_iter - 1)  + EVM_withDPD) / num_iter;

            CurrentRes.ACPR_noDPD.first = (CurrentRes.ACPR_noDPD.first * (num_iter - 1) + ACLR_noDPD.first) / num_iter;
            CurrentRes.ACPR_noDPD.second = (CurrentRes.ACPR_noDPD.second * (num_iter - 1) + ACLR_noDPD.second) / num_iter;

            CurrentRes.ACPR_withDPD.first = (CurrentRes.ACPR_withDPD.first * (num_iter - 1) + ACLR_withDPD.first) / num_iter;
            CurrentRes.ACPR_withDPD.second = (CurrentRes.ACPR_withDPD.second * (num_iter - 1) + ACLR_withDPD.second) / num_iter;

            CurrentRes.P_formed_noDPD = (CurrentRes.P_formed_noDPD * (num_iter - 1) + P_formed_noDPD) / num_iter;
            CurrentRes.P_formed_withDPD = (CurrentRes.P_formed_withDPD * (num_iter - 1)  + P_formed_withDPD) / num_iter;

            CurrentRes.P_emitted_noDPD = (CurrentRes.P_emitted_noDPD * (num_iter - 1) + P_emitted_noDPD) / num_iter;
            CurrentRes.P_emitted_withDPD = (CurrentRes.P_emitted_withDPD * (num_iter - 1)  + P_emitted_withDPD) / num_iter;

            CurrentRes.Gain_noDPD = (CurrentRes.Gain_noDPD * (num_iter - 1) + Gain_noDPD) / num_iter;
            CurrentRes.Gain_withDPD = (CurrentRes.Gain_withDPD * (num_iter - 1)  + Gain_withDPD) / num_iter;

            CurrentRes.PARP_noDPD = (CurrentRes.PARP_noDPD * (num_iter - 1)  + PARP_noDPD) / num_iter;
            CurrentRes.PARP_withDPD = (CurrentRes.PARP_withDPD * (num_iter - 1)  + PARP_withDPD) / num_iter;

            CurrentRes.NMSE_noDPD = (CurrentRes.NMSE_noDPD * (num_iter - 1)  + NMSE_noDPD) / num_iter;
            CurrentRes.NMSE_withDPD = (CurrentRes.NMSE_withDPD * (num_iter - 1)  + NMSE_withDPD) / num_iter;

            CurrentRes.CurrentOBO_noDPD = (CurrentRes.CurrentOBO_noDPD * (num_iter - 1)  + CurrentOBO_noDPD) / num_iter;
            CurrentRes.CurrentOBO_withDPD = (CurrentRes.CurrentOBO_withDPD * (num_iter - 1)  + CurrentOBO_withDPD) / num_iter;
        }
        else if(CurrentRecalcNeeds.CycleMode == false) {
            num_iter = 1;
            CurrentRes.ACPR_noDPD = MyMetricsEval.computeACPR(freq[0], PSDs[2], CurrentRes.BB, CurrentRes.BB * MySource.bb_delta, MySource.SigType);
            CurrentRes.ACPR_withDPD = MyMetricsEval.computeACPR(freq[0], PSDs[3], CurrentRes.BB, CurrentRes.BB * MySource.bb_delta, MySource.SigType);

            std::tie(CurrentRes.BER_noDPD, CurrentRes.BER_withDPD) = MyMetricsEval.Calc_BER(MySymbols);

            std::tie(CurrentRes.EVM_noDPD, CurrentRes.EVM_withDPD) = MyMetricsEval.Calc_EVM(MySymbols);

            CurrentRes.P_formed_noDPD = MyMetricsEval.compute_av_P_dBm(CurrentRes.tx_sig);
            CurrentRes.P_formed_withDPD = MyMetricsEval.compute_av_P_dBm(CurrentRes.tx_plus_dpd_sig);

            CurrentRes.P_emitted_noDPD = MyMetricsEval.compute_av_P_dBm(CurrentRes.pa_sig);
            CurrentRes.P_emitted_withDPD = MyMetricsEval.compute_av_P_dBm(CurrentRes.pa_plus_dpd_sig);

            CurrentRes.Gain_noDPD = CurrentRes.P_emitted_noDPD - CurrentRes.P_formed_noDPD;
            CurrentRes.Gain_withDPD = CurrentRes.P_emitted_withDPD - CurrentRes.P_formed_withDPD;

            CurrentRes.PARP_noDPD = MyMetricsEval.computePAPR_dB(CurrentRes.tx_sig);
            CurrentRes.PARP_withDPD = MyMetricsEval.computePAPR_dB(CurrentRes.tx_plus_dpd_sig);

            CurrentRes.NMSE_noDPD = MyMetricsEval.computeNMSE_dB(CurrentRes.tx_sig, CurrentRes.pa_sig);
            CurrentRes.NMSE_withDPD = MyMetricsEval.computeNMSE_dB(CurrentRes.tx_sig, CurrentRes.pa_plus_dpd_sig);

            CurrentRes.CurrentOBO_noDPD = MyMetricsEval.computeOBO_dB(CurrentRes.pa_sig, MySource.Pout_sat_dBm);
            CurrentRes.CurrentOBO_withDPD = MyMetricsEval.computeOBO_dB(CurrentRes.pa_plus_dpd_sig_noisy, MySource.Pout_sat_dBm);
        }
    }
    qDebug() << "Metrics time:" << timer.elapsed() << "ms";
}

void SignalProcessing::GeneratePacksOfSymbols(std::vector<Symbols>& Symbols, Source& source, NeedToRecalc& CurrentRecalcNeeds)
{
    if(CurrentRecalcNeeds.RecalcSymbols || CurrentRecalcNeeds.FullRecalc) {
        Symbols.clear();
        if(source.SigType == "FDMA")
            for(int i = 0; i < source.FDMA_num_subcarriers; i++) {
                Symbols.push_back(GenerateNSymbols(source));
                SymsAddNoise(source, Symbols[i].tr_sym_clean, Symbols[i].tr_sym_noisy);
            }
        else { Symbols.push_back(GenerateNSymbols(source)); SymsAddNoise(source, Symbols[0].tr_sym_clean, Symbols[0].tr_sym_noisy); }
    }
    else {
        if(source.SigType == "FDMA")
            for(int i = 0; i < source.FDMA_num_subcarriers; i++)
                SymsAddNoise(source, Symbols[i].tr_sym_clean, Symbols[i].tr_sym_noisy);
        else { SymsAddNoise(source, Symbols[0].tr_sym_clean, Symbols[0].tr_sym_noisy); }
    }
    CurrentRecalcNeeds.RecalcSig = true;
}

Symbols SignalProcessing::GenerateNSymbols(Source& source)
{
    Symbols temp;
    temp.resize(source.NumSym);

    int bitsPerSym = std::log2(source.M);
    int numBits = source.NumSym * bitsPerSym;

    temp.data_tx.resize(numBits);
    temp.data_rx.resize(numBits);
    temp.data_rx_with_DPD.resize(numBits);

    for (int i = 0; i < numBits; i++)
        temp.data_tx[i] = QRandomGenerator::global()->bounded(2); // 0 или 1

    for (int i = 0; i < source.NumSym; i++) {
        int idx = 0;

        for (int b = 0; b < bitsPerSym; b++) {
            idx = (idx << 1) | temp.data_tx[i * bitsPerSym + b];
        }

        temp.sym_idx_tx[i] = idx;
    }

    if(source.ModType == "BPSK") {
        for (int i = 0; i < source.NumSym; i++) {
            int s = temp.sym_idx_tx[i];
            temp.tr_sym_clean[i] = BPSK_const[s];
        }
    }
    else if(source.ModType == "QPSK") {
        for (int i = 0; i < source.NumSym; i++) {
            int s = temp.sym_idx_tx[i];
            temp.tr_sym_clean[i] = QPSK_const[s];
        }
    }
    else if(source.ModType == "16QAM") {
        for (int i = 0; i < source.NumSym; i++) {
            int s = temp.sym_idx_tx[i];
            temp.tr_sym_clean[i] = QAM16_const[s];
        }
    }
    else if(source.ModType == "64QAM") {
        for (int i = 0; i < source.NumSym; i++) {
            int s = temp.sym_idx_tx[i];
            temp.tr_sym_clean[i] = QAM64_const[s];
        }
    }
    else if(source.ModType == "APSK") {
        //int ring = ...;   // номер кольца
        //int idx  = ...;   // индекс на кольце
        //double r = RADi[ring];
        //double phi = 2 * M_PI * idx / Mring[ring];
        //
        //sym_clean[i] = std::polar(r, phi);
    }

    double P_signal = 0;
    for (auto& s : temp.tr_sym_clean)
        P_signal += std::norm(s);
    P_signal /= source.NumSym;

    double P_noise = P_signal / pow(10.0, source.SNRSymdB / 10.0);
    double sigma = std::sqrt(P_noise / 2);

    std::mt19937 rng(12345);
    std::normal_distribution<double> gauss(0.0, sigma);

    for (int i = 0; i < source.NumSym; i++) {
        temp.tr_sym_noisy[i] = temp.tr_sym_clean[i] +
                                     std::complex<double>(gauss(rng), gauss(rng));
    }

    return temp;
}

void SignalProcessing::SymsAddNoise(Source& source,
                                    std::vector<std::complex<double>>& symbols_clean, std::vector<std::complex<double>>& symbols_noisy) {
    double P_signal = 0;
    for (auto& s : symbols_clean)
        P_signal += std::norm(s);
    P_signal /= source.NumSym;

    double P_noise = P_signal / pow(10.0, source.SNRSymdB / 10.0);
    double sigma = std::sqrt(P_noise / 2);

    std::mt19937 rng(12345);
    std::normal_distribution<double> gauss(0.0, sigma);

    for (int i = 0; i < source.NumSym; i++) {
        symbols_noisy[i] = symbols_clean[i] +
                            std::complex<double>(gauss(rng), gauss(rng));
    }
}

int SignalProcessing::DemodulateSymbol(const std::complex<double> &r, const std::vector<std::complex<double> > &constellation)
{
    double minDist = std::numeric_limits<double>::max();
    int bestIndex = 0;

    for (size_t i = 0; i < constellation.size(); i++)
    {
        double dist = std::norm(r - constellation[i]);

        if (dist < minDist)
        {
            minDist = dist;
            bestIndex = i;
        }
    }

    return bestIndex;
}

void SignalProcessing::Demodulate(Symbols& symbols,const std::vector<std::complex<double>>& constellation, int M)
{
    int bitsPerSym = static_cast<int>(std::log2(M));
    int numBits = symbols.rec_sym_noisy.size() * bitsPerSym;

    symbols.data_rx.resize(numBits);
    symbols.data_rx_with_DPD.resize(numBits);

    for (size_t i = 0; i < symbols.rec_sym_noisy.size(); i++)
    {
        symbols.sym_idx_rx_noDPD[i] = DemodulateSymbol(symbols.rec_sym_noisy[i], constellation);
        symbols.sym_idx_rx_withDPD[i] = DemodulateSymbol(symbols.rec_sym_noisy_with_DPD[i], constellation);

        for (int b = 0; b < bitsPerSym; b++)
        {
            int shift = bitsPerSym - 1 - b;

            symbols.data_rx[i * bitsPerSym + b] =
                (symbols.sym_idx_rx_noDPD[i] >> shift) & 1;

            symbols.data_rx_with_DPD[i * bitsPerSym + b] =
                (symbols.sym_idx_rx_withDPD[i] >> shift) & 1;
        }
    }
}

void SignalProcessing::InitializeConstellations()
{
    BPSK_const =
        {
            {1,0},
            {-1,0}
        };

    double s = 1.0 / sqrt(2);
    QPSK_const =
        {
            { s,  s},
            {-s,  s},
            {-s, -s},
            { s, -s}
        };

    QAM16_const =
        {
            {-3,-3},{-3,-1},{-3,1},{-3,3},
            {-1,-3},{-1,-1},{-1,1},{-1,3},
            { 1,-3},{ 1,-1},{ 1,1},{ 1,3},
            { 3,-3},{ 3,-1},{ 3,1},{ 3,3}
        };
    double norm16 = 1.0 / std::sqrt(10.0);
    for(auto &s : QAM16_const)
        s *= norm16;

    QAM64_const =
        {
            {-7,-7},{-7,-5},{-7,-3},{-7,-1},{-7, 1},{-7, 3},{-7, 5},{-7, 7},
            {-5,-7},{-5,-5},{-5,-3},{-5,-1},{-5, 1},{-5, 3},{-5, 5},{-5, 7},
            {-3,-7},{-3,-5},{-3,-3},{-3,-1},{-3, 1},{-3, 3},{-3, 5},{-3, 7},
            {-1,-7},{-1,-5},{-1,-3},{-1,-1},{-1, 1},{-1, 3},{-1, 5},{-1, 7},

            { 1,-7},{ 1,-5},{ 1,-3},{ 1,-1},{ 1, 1},{ 1, 3},{ 1, 5},{ 1, 7},
            { 3,-7},{ 3,-5},{ 3,-3},{ 3,-1},{ 3, 1},{ 3, 3},{ 3, 5},{ 3, 7},
            { 5,-7},{ 5,-5},{ 5,-3},{ 5,-1},{ 5, 1},{ 5, 3},{ 5, 5},{ 5, 7},
            { 7,-7},{ 7,-5},{ 7,-3},{ 7,-1},{ 7, 1},{ 7, 3},{ 7, 5},{ 7, 7}
        };
    double norm64 = 1.0 / std::sqrt(42.0);
    for(auto &s : QAM64_const)
        s *= norm64;
}

void SignalProcessing::RecalcDPD(NeedToRecalc& CurrentRecalcNeeds)
{
    NeedToRecalc temp;
    temp.init();
    temp.RecalcNoiseSig = false;
    temp.RecalcNoiseSym = false;

    std::vector<Symbols> TrainSymbols;
    int sym = MySource.NumSym;
    MySource.NumSym = 2000;
    GlobalResults TrainRes;
    GeneratePacksOfSymbols(TrainSymbols, MySource, temp);
    TransmitSignalProcessing(MySource, TrainSymbols, temp, TrainRes);

    MyPAModels.ScaleToRMS_forPA(TrainRes.tx_sig, MySource);
    TrainRes.pa_sig = TrainRes.tx_sig;

    //if(MySource.IMUX_enabled)
    //    MyMux.apply(TrainRes.pa_sig, MuxKind::IMUX, TrainRes.BB, MySource.fs * MySource.oversampling);

    MyPAModels.ApplyPA(TrainRes.pa_sig, MySource);

    //if(MySource.OMUX_enabled)
    //    MyMux.apply(TrainRes.pa_sig, MuxKind::OMUX, TrainRes.BB, MySource.fs * MySource.oversampling);

    mydpd.train(TrainRes.tx_sig, TrainRes.pa_sig, MySource);
    MySource.NumSym = sym;
    CurrentRecalcNeeds.PARecalc = true;
}

std::vector<std::complex<double>>& SignalProcessing::getCurrentConstellation(Source& source)
{
    if(source.ModType == "BPSK") {
        return BPSK_const;
    }
    else if(source.ModType == "QPSK") {
        return QPSK_const;
    }
    else if(source.ModType == "16QAM") {
        return QAM16_const;
    }
    //else 64 QAM
    return QAM64_const;
}

OfdmParams SignalProcessing::GetOfdmParams(Source& source)
{
    OfdmParams ofdm_parms;
    ofdm_parms.CP = source.OFDM_cycle_prefix;
    ofdm_parms.GB_DC = source.OFDM_GB_DC;
    ofdm_parms.GB_Nyq = source.OFDM_GB_Nyq;
    ofdm_parms.Nfft = source.OFDM_Nfft;
    ofdm_parms.SNR_dB = source.SNRSig;
    ofdm_parms.fs = source.fs;
    ofdm_parms.oversampling = source.oversampling;
    return ofdm_parms;
}

ScParams SignalProcessing::GetSCParams(Source& source)
{
    ScParams sc_params;
    sc_params.SC_FilterType = source.SC_FilterType;
    sc_params.SC_filter_length = source.SC_filter_length;
    sc_params.SC_rolloff = source.SC_rolloff;
    sc_params.SC_symrate = source.SC_symrate;
    sc_params.SNR_dB = source.SNRSig;
    sc_params.fs = source.fs;
    sc_params.oversampling = source.oversampling;
    return sc_params;
}

FdmaParams SignalProcessing::GetFDMAParams(Source& source)
{
    FdmaParams fdma_params;
    fdma_params.FDMA_symrate = source.FDMA_symrate;
    fdma_params.FDMA_num_subcarriers = source.FDMA_num_subcarriers;
    fdma_params.FDMA_step_carrier = source.FDMA_step_carrier;
    fdma_params.SNRSig = source.SNRSig;
    fdma_params.fs = source.fs;
    fdma_params.oversampling = source.oversampling;
    fdma_params.EnableSNRSym = source.EnableSymSNR;

    return fdma_params;
}

void SignalProcessing::TransmitSignalProcessing(Source& source, std::vector<Symbols>& symbols, NeedToRecalc& CurrentRecalcNeeds, GlobalResults& CurRes)
{
    if(source.SigType == "OFDM") {
        OfdmParams ofdm_parms = GetOfdmParams(source);
        if(source.EnableSymSNR) CurrentOfdmResults = myOfdm.makeOfdm(symbols[0].tr_sym_noisy, ofdm_parms);
        else CurrentOfdmResults = myOfdm.makeOfdm(symbols[0].tr_sym_clean, ofdm_parms);
        CurRes.clear();
        CurRes.resize(CurrentOfdmResults.tx.size());
        CurRes.tx_sig = CurrentOfdmResults.tx;
        CurRes.time = CurrentOfdmResults.t;
        CurRes.BB = CurrentOfdmResults.BB;
    }
    else if(source.SigType == "FDMA") {
        FdmaParams fdma_parms = GetFDMAParams(source);
        ScParams sc_params = GetSCParams(source);
        CurrentFdmaResults = myFdma.generate(symbols, fdma_parms, sc_params);
        CurRes.clear();
        CurRes.resize(CurrentFdmaResults.tx.size());
        CurRes.tx_sig = CurrentFdmaResults.tx;
        CurRes.time = CurrentFdmaResults.t;
        CurRes.BB = CurrentFdmaResults.totalBandwidth;
    }
    else if(source.SigType == "SC") {
        ScParams sc_params = GetSCParams(source);
        if(source.EnableSymSNR) CurrentSCResults = mySC.makeSc(symbols[0].tr_sym_noisy, sc_params);
        else CurrentSCResults = mySC.makeSc(symbols[0].tr_sym_clean, sc_params);
        CurRes.clear();
        CurRes.resize(CurrentSCResults.tx.size());
        CurRes.tx_sig = CurrentSCResults.tx;
        CurRes.time = CurrentSCResults.t;
        CurRes.BB = CurrentSCResults.bandwidth;
    }
    CurrentRecalcNeeds.PARecalc = true;
}

void SignalProcessing::PAProcessing(Source& source, NeedToRecalc& CurrentRecalcNeeds, GlobalResults& CurRes)
{
    if(!CurRes.tx_sig.empty() && CurrentRecalcNeeds.PARecalc == true) {
        MyPAModels.ScaleToRMS_forPA(CurRes.tx_sig, source);
        CurRes.tx_plus_dpd_sig = CurRes.tx_sig;

        if(source.PredistorterType == "MP")
            CurRes.tx_plus_dpd_sig = mydpd.applyMP(CurRes.tx_plus_dpd_sig, source);
        else if(source.PredistorterType == "GMP")
            CurRes.tx_plus_dpd_sig = mydpd.applyGMP(CurRes.tx_plus_dpd_sig, source);

        //MyPAModels.ScaleToRMS_forPA(CurRes.tx_plus_dpd_sig, source);
        CurRes.pa_sig = CurRes.tx_sig;
        CurRes.pa_plus_dpd_sig = CurRes.tx_plus_dpd_sig;

        if(MySource.IMUX_enabled) {
            MyMux.apply(CurRes.pa_sig, MuxKind::IMUX, CurRes.BB, source.fs * source.oversampling);
            MyMux.apply(CurRes.pa_plus_dpd_sig, MuxKind::IMUX, CurRes.BB, source.fs * source.oversampling);
        }

        MyPAModels.ApplyPA(CurRes.pa_sig, MySource);
        MyPAModels.ApplyPA(CurRes.pa_plus_dpd_sig, MySource);

        if(MySource.OMUX_enabled) {
            MyMux.apply(CurRes.pa_sig, MuxKind::OMUX, CurRes.BB, source.fs * source.oversampling);
            MyMux.apply(CurRes.pa_plus_dpd_sig, MuxKind::OMUX, CurRes.BB, source.fs * source.oversampling);
        }

        CurRes.pa_sig_noisy = CurRes.pa_sig;
        CurRes.pa_plus_dpd_sig_noisy = CurRes.pa_plus_dpd_sig;
        addAwgn(CurRes.pa_sig_noisy, source.SNRSig);
        addAwgn(CurRes.pa_plus_dpd_sig_noisy, source.SNRSig);
        CurrentRecalcNeeds.RecRecalc = true;
    }
    else if(CurrentRecalcNeeds.RecalcNoiseSig) {
        CurRes.pa_sig_noisy = CurRes.pa_sig;
        CurRes.pa_plus_dpd_sig_noisy = CurRes.pa_plus_dpd_sig;
        addAwgn(CurRes.pa_sig_noisy, source.SNRSig);
        addAwgn(CurRes.pa_plus_dpd_sig_noisy, source.SNRSig);
        CurrentRecalcNeeds.RecRecalc = true;
    }
}

void SignalProcessing::ReceiveSignalProcessing(Source& source, std::vector<Symbols>& symbols, NeedToRecalc& CurrentRecalcNeeds, GlobalResults& CurRes)
{
    OfdmParams ofdm_parms = GetOfdmParams(source);
    ScParams sc_params = GetSCParams(source);
    FdmaParams fdma_params = GetFDMAParams(source);
    if(source.SigType == "OFDM") {
        symbols[0].rec_sym_noisy = myOfdm.ofdm_demodulate(CurRes.pa_sig_noisy, symbols[0].tr_sym_clean, ofdm_parms);
        symbols[0].rec_sym_noisy_with_DPD = myOfdm.ofdm_demodulate(CurRes.pa_plus_dpd_sig_noisy, symbols[0].tr_sym_clean, ofdm_parms);
    }
    else if(source.SigType == "SC") {
        symbols[0].rec_sym_noisy = mySC.demodulateSignal(CurRes.pa_sig_noisy, symbols[0].tr_sym_clean, sc_params, sc_params);
        symbols[0].rec_sym_noisy_with_DPD = mySC.demodulateSignal(CurRes.pa_plus_dpd_sig_noisy, symbols[0].tr_sym_clean, sc_params, sc_params);
    }
    else if(source.SigType == "FDMA") {
        std::vector<std::vector<std::complex<double>>> temp = myFdma.demodulate(CurRes.pa_sig_noisy, symbols, fdma_params, sc_params);
        std::vector<std::vector<std::complex<double>>> temp_with_DPD = myFdma.demodulate(CurRes.pa_plus_dpd_sig_noisy, symbols, fdma_params, sc_params);
        for(int i = 0; i < symbols.size(); ++i) {
            symbols[i].rec_sym_noisy = temp[i];
            symbols[i].rec_sym_noisy_with_DPD = temp_with_DPD[i];
        }
    }

    for(int i = 0; i < symbols.size(); ++i)
        Demodulate(symbols[i], getCurrentConstellation(source), source.M);

    CurrentRecalcNeeds.MetricsRecalc = true;
}

void SignalProcessing::CalcPaCurve()
{
    QString model_type = MySource.PAModel;

    if(model_type == "Wiener")
        model_type = MySource.W_StaticNonlinModel;
    else if(model_type == "Hammerstein")
        model_type = MySource.H_StaticNonlinModel;
    else if(model_type == "Wiener-Hammerstein")
        model_type = MySource.H_StaticNonlinModel;

    constexpr double eps = 1e-15;
    constexpr double rad2deg = 180.0 / M_PI;

    // -------------------------------------------------
    // Saturation points from GUI
    // -------------------------------------------------
    double Pin_sat_dBm  = MySource.Pin_sat_dBm;    // X saturation point
    double Pout_sat_dBm = MySource.Pout_sat_dBm;   // Y saturation point

    // -------------------------------------------------
    // AM/AM and AM/PM models
    // -------------------------------------------------
    auto calcAMAM = [&](double r) -> double
    {
        if(model_type == "Saleh")
        {
            return (MySource.SalehCoeffs[0] * r) /
                   (1.0 + MySource.SalehCoeffs[1] * qPow(r, 2.0));
        }
        else if(model_type == "Rapp")
        {
            return r /
                   qPow(1.0 + qPow(r / MySource.RappCoeffs[0],
                                   2.0 * MySource.RappCoeffs[1]),
                        1.0 / (2.0 * MySource.RappCoeffs[1]));
        }
        else if(model_type == "Ghorbani")
        {
            return MySource.GhorbaniCoeffs[0] * r /
                   (1.0
                    + MySource.GhorbaniCoeffs[1] * qPow(r, 2.0)
                    + MySource.GhorbaniCoeffs[2] * qPow(r, 4.0));
        }

        return r;
    };

    auto calcAMPM = [&](double r) -> double
    {
        if(model_type == "Saleh")
        {
            return (MySource.SalehCoeffs[2] * qPow(r, 2.0)) /
                   (1.0 + MySource.SalehCoeffs[3] * qPow(r, 2.0)) *
                   rad2deg;
        }
        else if(model_type == "Rapp")
        {
            return 0.0;
        }
        else if(model_type == "Ghorbani")
        {
            return MySource.GhorbaniCoeffs[3] * qPow(r, 2.0) /
                   (1.0
                    + MySource.GhorbaniCoeffs[4] * qPow(r, 2.0)
                    + MySource.GhorbaniCoeffs[5] * qPow(r, 4.0)) *
                   rad2deg;
        }

        return 0.0;
    };

    // -------------------------------------------------
    // 1. Raw AM/AM and AM/PM
    // -------------------------------------------------
    for(int i = 0; i < PACurve->point_num; i++)
    {
        double r = PACurve->r_in[i];

        PACurve->Aout[i] = calcAMAM(r);
        PACurve->Phi[i]  = calcAMPM(r);
    }

    // -------------------------------------------------
    // 2. Saturation point in the raw model
    // -------------------------------------------------
    auto maxItA = std::max_element(PACurve->Aout.begin(),
                                   PACurve->Aout.end());

    int maxIndex = std::distance(PACurve->Aout.begin(), maxItA);

    double Asat = std::max(*maxItA, eps);
    double r_sat = std::max(PACurve->r_in[maxIndex], eps);

    // -------------------------------------------------
    // 3. Absolute input power axis
    //
    // Pin_dBm = Pin_sat_dBm + 20log10(r / r_sat)
    // -------------------------------------------------
    for(int i = 0; i < PACurve->point_num; i++)
    {
        double r_norm = std::max(PACurve->r_in[i] / r_sat, eps);

        PACurve->P_in_abs.dB[i] =
            Pin_sat_dBm + 20.0 * std::log10(r_norm);

        PACurve->P_in_abs.linear[i] =
            qPow(10.0, PACurve->P_in_abs.dB[i] / 10.0);
    }

    // -------------------------------------------------
    // 4. Absolute output power axis
    //
    // Pout_dBm = Pout_sat_dBm + 20log10(A / Asat)
    // -------------------------------------------------
    for(int i = 0; i < PACurve->point_num; i++)
    {
        double A_norm = std::max(PACurve->Aout[i] / Asat, eps);

        PACurve->P_out_abs.dB[i] =
            Pout_sat_dBm + 20.0 * std::log10(A_norm);

        PACurve->P_out_abs.linear[i] =
            qPow(10.0, PACurve->P_out_abs.dB[i] / 10.0);
    }

    // -------------------------------------------------
    // 5. Normalized input/output curves
    // -------------------------------------------------
    for(int i = 0; i < PACurve->point_num; i++)
    {
        PACurve->P_in_norm.dB[i] =
            PACurve->P_in_abs.dB[i] - Pin_sat_dBm;

        PACurve->P_in_norm.linear[i] =
            qPow(10.0, PACurve->P_in_norm.dB[i] / 10.0);

        PACurve->P_out_norm.dB[i] =
            PACurve->P_out_abs.dB[i] - Pout_sat_dBm;

        PACurve->P_out_norm.linear[i] =
            qPow(10.0, PACurve->P_out_norm.dB[i] / 10.0);
    }

    // -------------------------------------------------
    // 6. Working point from IBO
    //
    // IBO = Pin_sat / Pin_work
    // Pin_work_dBm = Pin_sat_dBm - IBO_dB
    // r_work = r_sat / sqrt(IBO_linear)
    // -------------------------------------------------
    double IBO_linear = qPow(10.0, MySource.IBO_dB / 10.0);
    double r_work = r_sat / std::sqrt(std::max(IBO_linear, eps));

    double A_work = calcAMAM(r_work);
    double Phi_work = calcAMPM(r_work);

    double A_work_norm = std::max(A_work / Asat, eps);

    PACurve->Phi_work_grad[0] = Phi_work;

    // -------------------------------------------------
    // 7. Absolute working point
    // -------------------------------------------------
    PACurve->Working_point_dB_abs.x[0] =
        Pin_sat_dBm - MySource.IBO_dB;

    PACurve->Working_point_dB_abs.y[0] =
        Pout_sat_dBm + 20.0 * std::log10(A_work_norm);

    PACurve->Working_point_linear_abs.x[0] =
        qPow(10.0, PACurve->Working_point_dB_abs.x[0] / 10.0);

    PACurve->Working_point_linear_abs.y[0] =
        qPow(10.0, PACurve->Working_point_dB_abs.y[0] / 10.0);

    // -------------------------------------------------
    // 8. Normalized working point
    // -------------------------------------------------
    PACurve->Working_point_dB_norm.x[0] =
        -MySource.IBO_dB;

    PACurve->Working_point_dB_norm.y[0] =
        20.0 * std::log10(A_work_norm);

    PACurve->Working_point_linear_norm.x[0] =
        qPow(10.0, PACurve->Working_point_dB_norm.x[0] / 10.0);

    PACurve->Working_point_linear_norm.y[0] =
        qPow(10.0, PACurve->Working_point_dB_norm.y[0] / 10.0);
}

void SignalProcessing::ApplyFIRWithMemory(std::vector<double>& amplitude, std::vector<double>& phase,
    const std::vector<double>& FIR_Coefs, int numTaps)
{
    const size_t N = amplitude.size();

    double C = FIR_Coefs[0];
    double alpha = FIR_Coefs[1];

    // Формируем коэффициенты FIR
    std::vector<double> h(numTaps);
    for (int m = 0; m < numTaps; ++m)
    {
        h[m] = C * std::pow(alpha, m);
    }

    // Нормировка (чтобы не менять общий gain)
    double sum = 0.0;
    for (double val : h)
        sum += val;

    for (double& val : h)
        val /= sum;

    // Формируем комплексный вход
    std::vector<std::complex<double>> input(N);
    for (size_t n = 0; n < N; ++n)
    {
        input[n] = std::polar(amplitude[n], phase[n]);
    }

    // FIR
    std::vector<std::complex<double>> output(N, {0.0, 0.0});

    for (size_t n = 0; n < N; ++n)
    {
        for (int m = 0; m < numTaps; ++m)
        {
            if (n >= static_cast<size_t>(m))
            {
                output[n] += h[m] * input[n - m];
            }
        }
    }

    // Обратно в амплитуду и фазу
    for (size_t n = 0; n < N; ++n)
    {
        amplitude[n] = std::abs(output[n]);
        phase[n] = 0; //std::arg(output[n]);
    }
}

void SignalProcessing::DataUpdate(Source &UISource)
{
    MySource.ModType = UISource.ModType;
    MySource.SNRSymdB = UISource.SNRSymdB;
    MySource.NumSym = UISource.NumSym;
    MySource.M = UISource.M;
    MySource.EnableSymSNR = UISource.EnableSymSNR;

    MySource.SigType = UISource.SigType;
    MySource.SC_symrate = UISource.SC_symrate;
    MySource.SC_rolloff = UISource.SC_rolloff;
    MySource.SC_filter_length = UISource.SC_filter_length;
    MySource.SC_FilterType = UISource.SC_FilterType;

    MySource.OFDM_Nfft = UISource.OFDM_Nfft;
    MySource.OFDM_Nfft = UISource.OFDM_Nfft;
    MySource.OFDM_GB_DC = UISource.OFDM_GB_DC;
    MySource.OFDM_GB_Nyq = UISource.OFDM_GB_Nyq;
    MySource.OFDM_cycle_prefix = UISource.OFDM_cycle_prefix;

    MySource.FDMA_symrate = UISource.FDMA_symrate;
    MySource.FDMA_num_subcarriers = UISource.FDMA_num_subcarriers;
    MySource.FDMA_step_carrier = UISource.FDMA_step_carrier;

    MySource.oversampling = UISource.oversampling;
    MySource.fs = UISource.fs;
    MySource.SNRSig = UISource.SNRSig;
    MySource.bb_delta = UISource.bb_delta;

    MySource.Pin_sat_dBm = UISource.Pin_sat_dBm;
    MySource.Pout_sat_dBm = UISource.Pout_sat_dBm;
    MySource.PAModel = UISource.PAModel;
    MySource.IBO_dB = UISource.IBO_dB;
    MySource.linear_gain_dB = UISource.linear_gain_dB;
    MySource.SalehCoeffs = UISource.SalehCoeffs;
    MySource.RappCoeffs = UISource.RappCoeffs;
    MySource.GhorbaniCoeffs = UISource.GhorbaniCoeffs;
    MySource.W_FIRCoeffs = UISource.W_FIRCoeffs;
    MySource.W_StaticNonlinModel = UISource.W_StaticNonlinModel;
    MySource.H_FIRCoeffs = UISource.H_FIRCoeffs;
    MySource.H_StaticNonlinModel = UISource.H_StaticNonlinModel;
    MySource.WH_FIRCoeffs = UISource.WH_FIRCoeffs;
    MySource.WH_StaticNonlinModel = UISource.WH_StaticNonlinModel;

    MySource.MP_M = UISource.MP_M;
    MySource.MP_P = UISource.MP_P;

    MySource.GMP_M = UISource.GMP_M;
    MySource.GMP_P = UISource.GMP_P;
    MySource.GMP_L_lag = UISource.GMP_L_lag;
    MySource.GMP_L_lead = UISource.GMP_L_lead;

    MySource.NormalizationType = UISource.NormalizationType;
    MySource.PredistorterType = UISource.PredistorterType;
    MySource.DPDAutoRecalc = UISource.DPDAutoRecalc;
    MySource.Enable_even_P = UISource.Enable_even_P;

    MySource.IMUX_enabled = UISource.IMUX_enabled;
    MySource.OMUX_enabled = UISource.OMUX_enabled;
}

void SignalProcessing::addAwgn(std::vector<std::complex<double>> &x, double SNR_dB)
{
    // Считаем мощность сигнала
    double power = 0;
    for(const auto& v : x)
        power += std::norm(v);  // norm = real^2 + imag^2
    power /= x.size();

    // Считаем дисперсию шума
    double snr = std::pow(10.0, SNR_dB / 10.0);
    double noiseVar = power / snr;
    double noiseStd = std::sqrt(noiseVar / 2);  // /2 потому что шум комплексный (I и Q компоненты)

    // Генератор шума
    static std::default_random_engine gen(std::random_device{}());  // static чтоб каждый раз не создавать
    std::normal_distribution<double> dist(0.0, noiseStd);

    // Добавляем шум
    for(int i = 0; i < x.size(); ++i) {
        double noiseI = dist(gen);  // действительная часть
        double noiseQ = dist(gen);  // мнимая часть
        x[i] += std::complex<double>(noiseI, noiseQ);
    }
}

Symbols &SignalProcessing::getSymbols()
{
        return MySymbols[0];
}

OfdmResult SignalProcessing::getOfdmRes()
{
    return CurrentOfdmResults;
}

ScResult SignalProcessing::getSCRes()
{
    return CurrentSCResults;
}

FdmaResult SignalProcessing::getFDMARes()
{
    return CurrentFdmaResults;
}

GlobalResults &SignalProcessing::getTimeSignal()
{
        return CurrentRes;
}

std::vector<std::vector<double>>& SignalProcessing::getFreq()
{
    return freq;
}

std::vector<std::vector<double>>& SignalProcessing::getPSDs()
{
    return PSDs;
}

ImuxOmux &SignalProcessing::getMux()
{
    return MyMux;
}

PaCurve &SignalProcessing::getPaCurve()
{
    return *PACurve;
}

void SignalProcessing::clear_OFDM_buffs()
{
    OFDM_tx_sig_buff.clear();
    OFDM_tx_plus_DPD_sig_buff.clear();
    OFDM_pa_sig_buff.clear();
    OFDM_pa_plus_DPD_sig_buff.clear();
}

