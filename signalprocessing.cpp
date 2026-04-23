#include "signalprocessing.h"

SignalProcessing::SignalProcessing() : myFdma(mySC), mydpd() {
    PACurve = new PaCurve(200);
    // Input power in dB and linear
    double dB_start = -25;
    double dB_end = 5;
    double dB_step = (dB_end - dB_start) / (PACurve->point_num - 1);
    for(int i = 0; i < PACurve->point_num; i ++) {
        PACurve->P_in_abs.dB[i] = dB_start + i * dB_step;
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
            P_formed_noDPD = MyMetricsEval.compute_av_P(CurrentRes.tx_sig);
            P_formed_withDPD = MyMetricsEval.compute_av_P(CurrentRes.tx_plus_dpd_sig);
            P_emitted_noDPD = MyMetricsEval.compute_av_P(CurrentRes.pa_sig);
            P_emitted_withDPD = MyMetricsEval.compute_av_P(CurrentRes.pa_plus_dpd_sig);
            Gain_noDPD = MyMetricsEval.compute_av_P_G(CurrentRes.P_formed_noDPD, CurrentRes.P_emitted_noDPD);
            Gain_withDPD = MyMetricsEval.compute_av_P_G(CurrentRes.P_formed_noDPD, CurrentRes.P_emitted_withDPD);
            PARP_noDPD = MyMetricsEval.computePAPR_dB(CurrentRes.tx_sig);
            PARP_withDPD = MyMetricsEval.computePAPR_dB(CurrentRes.tx_plus_dpd_sig);
            NMSE_noDPD = MyMetricsEval.computeNMSE_dB(CurrentRes.tx_sig, CurrentRes.pa_sig);
            NMSE_withDPD = MyMetricsEval.computeNMSE_dB(CurrentRes.tx_sig, CurrentRes.pa_plus_dpd_sig);
            CurrentOBO_noDPD = 10 * std::log10(MySource.PA_sat / CurrentRes.P_emitted_noDPD);
            CurrentOBO_withDPD = 10 * std::log10(MySource.PA_sat / CurrentRes.P_emitted_withDPD);

            ACLR_noDPD = MyMetricsEval.computeACPR(freq[0], PSDs[1], CurrentRes.BB, CurrentRes.BB * MySource.bb_delta, MySource.SigType);
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
            CurrentRes.ACPR_noDPD = MyMetricsEval.computeACPR(freq[0], PSDs[1], CurrentRes.BB, CurrentRes.BB * MySource.bb_delta, MySource.SigType);
            CurrentRes.ACPR_withDPD = MyMetricsEval.computeACPR(freq[0], PSDs[3], CurrentRes.BB, CurrentRes.BB * MySource.bb_delta, MySource.SigType);

            std::tie(CurrentRes.BER_noDPD, CurrentRes.BER_withDPD) = MyMetricsEval.Calc_BER(MySymbols);

            std::tie(CurrentRes.EVM_noDPD, CurrentRes.EVM_withDPD) = MyMetricsEval.Calc_EVM(MySymbols);

            CurrentRes.P_formed_noDPD = MyMetricsEval.compute_av_P(CurrentRes.tx_sig);
            CurrentRes.P_formed_withDPD = MyMetricsEval.compute_av_P(CurrentRes.tx_plus_dpd_sig);

            CurrentRes.P_emitted_noDPD = MyMetricsEval.compute_av_P(CurrentRes.pa_sig);
            CurrentRes.P_emitted_withDPD = MyMetricsEval.compute_av_P(CurrentRes.pa_plus_dpd_sig);

            CurrentRes.Gain_noDPD = MyMetricsEval.compute_av_P_G(CurrentRes.P_formed_noDPD, CurrentRes.P_emitted_noDPD);
            CurrentRes.Gain_withDPD = MyMetricsEval.compute_av_P_G(CurrentRes.P_formed_noDPD, CurrentRes.P_emitted_withDPD);

            CurrentRes.PARP_noDPD = MyMetricsEval.computePAPR_dB(CurrentRes.tx_sig);
            CurrentRes.PARP_withDPD = MyMetricsEval.computePAPR_dB(CurrentRes.tx_plus_dpd_sig);

            CurrentRes.NMSE_noDPD = MyMetricsEval.computeNMSE_dB(CurrentRes.tx_sig, CurrentRes.pa_sig);
            CurrentRes.NMSE_withDPD = MyMetricsEval.computeNMSE_dB(CurrentRes.tx_sig, CurrentRes.pa_plus_dpd_sig);

            CurrentRes.CurrentOBO_noDPD = 10 * std::log10(MySource.PA_sat / CurrentRes.P_emitted_noDPD);
            CurrentRes.CurrentOBO_withDPD = 10 * std::log10(MySource.PA_sat / CurrentRes.P_emitted_withDPD);
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
    temp.tr_sym_clean.resize(source.NumSym);
    temp.tr_sym_noisy.resize(source.NumSym);
    temp.data_tx.resize(source.NumSym);
    for (int i = 0; i < source.NumSym; i++)
        temp.data_tx[i] = QRandomGenerator::global()->bounded(source.M);

    if(source.ModType == "BPSK") {
        for (int i = 0; i < source.NumSym; i++) {
            int s = temp.data_tx[i];
            temp.tr_sym_clean[i] = BPSK_const[s];
        }
    }
    else if(source.ModType == "QPSK") {
        for (int i = 0; i < source.NumSym; i++) {
            int s = temp.data_tx[i];
            temp.tr_sym_clean[i] = QPSK_const[s];
        }
    }
    else if(source.ModType == "16QAM") {
        for (int i = 0; i < source.NumSym; i++) {
            int s = temp.data_tx[i];
            temp.tr_sym_clean[i] = QAM16_const[s];
        }
    }
    else if(source.ModType == "64QAM") {
        for (int i = 0; i < source.NumSym; i++) {
            int s = temp.data_tx[i];
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

void SignalProcessing::Demodulate(Symbols& symbols, const std::vector<std::complex<double> > &constellation)
{
    symbols.data_rx.resize(symbols.data_tx.size());
    symbols.data_rx_with_DPD.resize(symbols.data_tx.size());
    for (size_t i = 0; i < symbols.rec_sym_noisy.size(); i++)
    {
        symbols.data_rx[i] = DemodulateSymbol(symbols.rec_sym_noisy[i], constellation);
        symbols.data_rx_with_DPD[i] = DemodulateSymbol(symbols.rec_sym_noisy_with_DPD[i], constellation);
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

    if(MySource.PAModel == "Saleh")
        MyPAModels.SalehModel(TrainRes.pa_sig, MySource.SalehCoeffs, MySource.linear_gain_dB, MySource.IBO_dB);
    else if (MySource.PAModel == "Rapp")
        MyPAModels.RappModel(TrainRes.pa_sig, MySource.RappCoeffs, MySource.linear_gain_dB, MySource.IBO_dB);
    else if (MySource.PAModel == "Ghorbani")
        MyPAModels.GhorbaniModel(TrainRes.pa_sig, MySource.GhorbaniCoeffs, MySource.linear_gain_dB, MySource.IBO_dB);
    else if (MySource.PAModel == "Wiener") {
        if(MySource.W_StaticNonlinModel == "Saleh")
            MyPAModels.WienerModel(TrainRes.pa_sig, MySource.W_StaticNonlinModel,
                                   MySource.SalehCoeffs, MySource.W_FIRCoeffs, MySource.linear_gain_dB, MySource.IBO_dB);
        else if(MySource.W_StaticNonlinModel == "Rapp")
            MyPAModels.WienerModel(TrainRes.pa_sig, MySource.W_StaticNonlinModel,
                                   MySource.RappCoeffs, MySource.W_FIRCoeffs, MySource.linear_gain_dB, MySource.IBO_dB);
        else if(MySource.W_StaticNonlinModel == "Ghorbani")
            MyPAModels.WienerModel(TrainRes.pa_sig, MySource.W_StaticNonlinModel,
                                   MySource.GhorbaniCoeffs, MySource.W_FIRCoeffs, MySource.linear_gain_dB, MySource.IBO_dB);
    }
    else if (MySource.PAModel == "Hammerstein") {
        if(MySource.H_StaticNonlinModel == "Saleh")
            MyPAModels.HammersteinModel(TrainRes.pa_sig, MySource.H_StaticNonlinModel,
                                   MySource.SalehCoeffs, MySource.H_FIRCoeffs, MySource.linear_gain_dB, MySource.IBO_dB);
        else if(MySource.H_StaticNonlinModel == "Rapp")
            MyPAModels.HammersteinModel(TrainRes.pa_sig, MySource.H_StaticNonlinModel,
                                   MySource.RappCoeffs, MySource.H_FIRCoeffs, MySource.linear_gain_dB, MySource.IBO_dB);
        else if(MySource.H_StaticNonlinModel == "Ghorbani")
            MyPAModels.HammersteinModel(TrainRes.pa_sig, MySource.H_StaticNonlinModel,
                                   MySource.GhorbaniCoeffs, MySource.H_FIRCoeffs, MySource.linear_gain_dB, MySource.IBO_dB);
    }
    else if (MySource.PAModel == "Wiener-Hammerstein") {
        if(MySource.WH_StaticNonlinModel == "Saleh")
            MyPAModels.WHModel(TrainRes.pa_sig, MySource.WH_StaticNonlinModel,
                                   MySource.SalehCoeffs, MySource.WH_FIRCoeffs, MySource.linear_gain_dB, MySource.IBO_dB);
        else if(MySource.WH_StaticNonlinModel == "Rapp")
            MyPAModels.WHModel(TrainRes.pa_sig, MySource.WH_StaticNonlinModel,
                                   MySource.RappCoeffs, MySource.WH_FIRCoeffs, MySource.linear_gain_dB, MySource.IBO_dB);
        else if(MySource.WH_StaticNonlinModel == "Ghorbani")
            MyPAModels.WHModel(TrainRes.pa_sig, MySource.WH_StaticNonlinModel,
                                   MySource.GhorbaniCoeffs, MySource.WH_FIRCoeffs, MySource.linear_gain_dB, MySource.IBO_dB);
    }

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
    ofdm_parms.fc = source.OFDM_f_carrier;
    ofdm_parms.fs = source.fs;
    ofdm_parms.oversampling = source.oversampling;
    return ofdm_parms;
}

ScParams SignalProcessing::GetSCParams(Source& source)
{
    ScParams sc_params;
    sc_params.SC_FilterType = source.SC_FilterType;
    sc_params.SC_filter_length = source.SC_filter_length;
    sc_params.SC_f_carrier = source.SC_f_carrier;
    sc_params.SC_rolloff = source.SC_rolloff;
    sc_params.SC_symrate = source.SC_symrate;
    sc_params.SNR_dB = source.SNRSig;
    sc_params.fs = source.fs;
    sc_params.fc = source.SC_f_carrier;
    sc_params.oversampling = source.oversampling;
    return sc_params;
}

FdmaParams SignalProcessing::GetFDMAParams(Source& source)
{
    FdmaParams fdma_params;
    fdma_params.FDMA_f_carrier = source.FDMA_f_carrier;
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
        CurRes.pa_sig = CurRes.tx_sig;
        CurRes.tx_plus_dpd_sig = CurRes.tx_sig;

        if(source.PredistorterType == "MP")
            CurRes.tx_plus_dpd_sig = mydpd.applyMP(CurRes.tx_plus_dpd_sig, source);
        else if(source.PredistorterType == "GMP")
            CurRes.tx_plus_dpd_sig = mydpd.applyGMP(CurRes.tx_plus_dpd_sig, source);

        MyPAModels.ScaleToRMS_forPA(CurRes.tx_plus_dpd_sig, source);
        CurRes.pa_plus_dpd_sig = CurRes.tx_plus_dpd_sig;

        if(source.PAModel == "Saleh") {
            MyPAModels.SalehModel(CurRes.pa_sig, source.SalehCoeffs, source.linear_gain_dB, source.IBO_dB);
            MyPAModels.SalehModel(CurRes.pa_plus_dpd_sig, source.SalehCoeffs, source.linear_gain_dB, source.IBO_dB);
        }
        else if (source.PAModel == "Rapp") {
            MyPAModels.RappModel(CurRes.pa_sig, source.RappCoeffs, source.linear_gain_dB, source.IBO_dB);
            MyPAModels.RappModel(CurRes.pa_plus_dpd_sig, source.RappCoeffs, source.linear_gain_dB, source.IBO_dB);
        }
        else if (source.PAModel == "Ghorbani") {
            MyPAModels.GhorbaniModel(CurRes.pa_sig, source.GhorbaniCoeffs, source.linear_gain_dB, source.IBO_dB);
            MyPAModels.GhorbaniModel(CurRes.pa_plus_dpd_sig, source.GhorbaniCoeffs, source.linear_gain_dB, source.IBO_dB);
        }
        else if (source.PAModel == "Wiener") {
            if(source.W_StaticNonlinModel == "Saleh") {
                MyPAModels.WienerModel(CurRes.pa_sig, source.W_StaticNonlinModel,
                                       source.SalehCoeffs, source.W_FIRCoeffs, source.linear_gain_dB, source.IBO_dB);
                MyPAModels.WienerModel(CurRes.pa_plus_dpd_sig, source.W_StaticNonlinModel,
                                       source.SalehCoeffs, source.W_FIRCoeffs, source.linear_gain_dB, source.IBO_dB);
            }
            else if(source.W_StaticNonlinModel == "Rapp") {
                MyPAModels.WienerModel(CurRes.pa_sig, source.W_StaticNonlinModel,
                                       source.RappCoeffs, source.W_FIRCoeffs, source.linear_gain_dB, source.IBO_dB);
                MyPAModels.WienerModel(CurRes.pa_plus_dpd_sig, source.W_StaticNonlinModel,
                                       source.RappCoeffs, source.W_FIRCoeffs, source.linear_gain_dB, source.IBO_dB);
            }
            else if(source.W_StaticNonlinModel == "Ghorbani") {
                MyPAModels.WienerModel(CurRes.pa_sig, source.W_StaticNonlinModel,
                                       source.GhorbaniCoeffs, source.W_FIRCoeffs, source.linear_gain_dB, source.IBO_dB);
                MyPAModels.WienerModel(CurRes.pa_plus_dpd_sig, source.W_StaticNonlinModel,
                                       source.GhorbaniCoeffs, source.W_FIRCoeffs, source.linear_gain_dB, source.IBO_dB);
            }
        }
        else if (source.PAModel == "Hammerstein") {
            if(source.H_StaticNonlinModel == "Saleh") {
                MyPAModels.HammersteinModel(CurRes.pa_sig, source.H_StaticNonlinModel,
                                       source.SalehCoeffs, source.H_FIRCoeffs, source.linear_gain_dB, source.IBO_dB);
                MyPAModels.HammersteinModel(CurRes.pa_plus_dpd_sig, source.H_StaticNonlinModel,
                                       source.SalehCoeffs, source.H_FIRCoeffs, source.linear_gain_dB, source.IBO_dB);
            }
            else if(source.H_StaticNonlinModel == "Rapp") {
                MyPAModels.HammersteinModel(CurRes.pa_sig, source.H_StaticNonlinModel,
                                       source.RappCoeffs, source.H_FIRCoeffs, source.linear_gain_dB, source.IBO_dB);
                MyPAModels.HammersteinModel(CurRes.pa_plus_dpd_sig, source.H_StaticNonlinModel,
                                       source.RappCoeffs, source.H_FIRCoeffs, source.linear_gain_dB, source.IBO_dB);
            }
            else if(source.H_StaticNonlinModel == "Ghorbani") {
                MyPAModels.HammersteinModel(CurRes.pa_sig, source.H_StaticNonlinModel,
                                       source.GhorbaniCoeffs, source.H_FIRCoeffs, source.linear_gain_dB, source.IBO_dB);
                MyPAModels.HammersteinModel(CurRes.pa_plus_dpd_sig, source.H_StaticNonlinModel,
                                       source.GhorbaniCoeffs, source.H_FIRCoeffs, source.linear_gain_dB, source.IBO_dB);
            }
        }
        else if (source.PAModel == "Wiener-Hammerstein") {
            if(source.WH_StaticNonlinModel == "Saleh") {
                MyPAModels.WHModel(CurRes.pa_sig, source.WH_StaticNonlinModel,
                                       source.SalehCoeffs, source.WH_FIRCoeffs, source.linear_gain_dB, source.IBO_dB);
                MyPAModels.WHModel(CurRes.pa_plus_dpd_sig, source.WH_StaticNonlinModel,
                                       source.SalehCoeffs, source.WH_FIRCoeffs, source.linear_gain_dB, source.IBO_dB);
            }
            else if(source.WH_StaticNonlinModel == "Rapp") {
                MyPAModels.WHModel(CurRes.pa_sig, source.WH_StaticNonlinModel,
                                       source.RappCoeffs, source.WH_FIRCoeffs, source.linear_gain_dB, source.IBO_dB);
                MyPAModels.WHModel(CurRes.pa_plus_dpd_sig, source.WH_StaticNonlinModel,
                                       source.RappCoeffs, source.WH_FIRCoeffs, source.linear_gain_dB, source.IBO_dB);
            }
            else if(source.WH_StaticNonlinModel == "Ghorbani") {
                MyPAModels.WHModel(CurRes.pa_sig, source.WH_StaticNonlinModel,
                                       source.GhorbaniCoeffs, source.WH_FIRCoeffs, source.linear_gain_dB, source.IBO_dB);
                MyPAModels.WHModel(CurRes.pa_plus_dpd_sig, source.WH_StaticNonlinModel,
                                       source.GhorbaniCoeffs, source.WH_FIRCoeffs, source.linear_gain_dB, source.IBO_dB);
            }
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
        Demodulate(symbols[i], getCurrentConstellation(source));

    CurrentRecalcNeeds.MetricsRecalc = true;
}

void SignalProcessing::CalcPaCurve()
{
    bool FIR_enable;
    QString model_type = MySource.PAModel;
    if(model_type == "Wiener") {
        model_type = MySource.W_StaticNonlinModel;
        FIR_enable = true;
    }

    double linear_gain = qPow(10.0, MySource.linear_gain_dB / 20);
    if(model_type == "Saleh") {
        for(int i = 0; i < PACurve->point_num; i ++) {
            // AM/AM (unnormalized)
            PACurve->Aout[i] = (MySource.SalehCoeffs[0] * PACurve->r_in[i]) / (1 + MySource.SalehCoeffs[1] * qPow(PACurve->r_in[i], 2.0)) ;
            // AM/PM
            PACurve->Phi[i] = (MySource.SalehCoeffs[2] * qPow(PACurve->r_in[i], 2)) /  (1 + MySource.SalehCoeffs[3] * qPow(PACurve->r_in[i], 2)) * 180 / 3.14;
        }
    }
    else if (model_type == "Rapp") {
        for(int i = 0; i < PACurve->point_num; i ++) {
            // AM/AM (unnormalized)
            PACurve->Aout[i] = linear_gain * PACurve->r_in[i] / qPow(1 + qPow(PACurve->r_in[i] / MySource.RappCoeffs[0], 2 * MySource.RappCoeffs[1]), 1 / (2 * MySource.RappCoeffs[1]));
            // AM/PM
            PACurve->Phi[i] = 0;
        }
    }
    else if (model_type == "Ghorbani") {
        for(int i = 0; i < PACurve->point_num; i ++) {
            // AM/AM (unnormalized)
            PACurve->Aout[i] = MySource.GhorbaniCoeffs[0] * PACurve->r_in[i] / (1 + MySource.GhorbaniCoeffs[1] * qPow(PACurve->r_in[i], 2) + MySource.GhorbaniCoeffs[2] * qPow(PACurve->r_in[i], 4));
            // AM/PM
            PACurve->Phi[i] = MySource.GhorbaniCoeffs[3] * qPow(PACurve->r_in[i], 2) / (1 + MySource.GhorbaniCoeffs[4] * qPow(PACurve->r_in[i], 2) + MySource.GhorbaniCoeffs[5] * qPow(PACurve->r_in[i], 4)) * 180 / 3.14;
        }
    }

    // Output power in dB and linear
    for(int i = 0; i < PACurve->point_num; i ++) {
        // P out in dB with linear gain
        PACurve->P_out_abs.dB[i] =  10 * std::log10(qPow(PACurve->Aout[i], 2)) + MySource.linear_gain_dB;
        PACurve->P_out_abs.linear[i] = qPow(10.0, PACurve->P_out_abs.dB[i] / 10.0);
    }

    // Make artificial saturation
    // --- ЭМПИРИЧЕСКАЯ сатурация (по факту, без аналитической асимптоты) ---
    auto maxIt = std::max_element(PACurve->P_out_abs.dB.begin(), PACurve->P_out_abs.dB.end());
    double Psat = *maxIt;               // значение максимума
    int maxIndex = maxIt - PACurve->P_out_abs.dB.begin();    // индекс максимума
    double Asat = PACurve->Aout[maxIndex];
    double P_in_sat_dB = PACurve->P_in_abs.dB[maxIndex];
    double P_in_sat_linear = PACurve->P_in_abs.linear[maxIndex];
    double r_sat = PACurve->r_in[maxIndex];

    //Нормированная выходная мощность: вычитаем Psat (gain_linear автоматически уйдёт)
    for (int i = 0; i < PACurve->point_num; i++) {
        PACurve->P_out_norm.dB[i] = PACurve->P_out_abs.dB[i] - Psat;
        PACurve->P_out_norm.linear[i] = qPow(10.0, PACurve->P_out_norm.dB[i] / 10.0);
    }

    // Нормированный вход: Pin/Pin_sat (в дБ и линейно) — используем эмпирическое r_sat
    // Замечание: т.к. r_in — амплитуда, используем 20*log10 для dB амплитуды
    for(int i = 0; i < PACurve->point_num; i++) {
        PACurve->P_in_norm.dB[i] = PACurve->P_in_abs.dB[i] - P_in_sat_dB;
        PACurve->P_in_norm.linear[i] = PACurve->P_in_abs.linear[i] / P_in_sat_linear;
    }

    // Рабочая точка
    double IBO_linear = qPow(10.0, MySource.IBO_dB / 10.0);
    double r_work = r_sat / std::sqrt(IBO_linear);

    double A_work;
    if(model_type == "Saleh") {
        A_work = (MySource.SalehCoeffs[0] * r_work) / (1 + MySource.SalehCoeffs[1] * qPow(r_work, 2));
        PACurve->Phi_work_grad[0] = (MySource.SalehCoeffs[2] * qPow(r_work, 2)) / (1 + MySource.SalehCoeffs[3] * qPow(r_work, 2)) * 180 / 3.14;
    }
    else if(model_type == "Rapp") {
        A_work = r_work / qPow(1 + qPow(r_work / MySource.RappCoeffs[0], 2 * MySource.RappCoeffs[1]), 1 / (2 * MySource.RappCoeffs[1]));
        PACurve->Phi_work_grad[0] = 0;
    }
    else if(model_type == "Ghorbani") {
        A_work = MySource.GhorbaniCoeffs[0] * r_work / (1 + MySource.GhorbaniCoeffs[1] * qPow(r_work, 2) + MySource.GhorbaniCoeffs[2] * qPow(r_work, 4));
        PACurve->Phi_work_grad[0] = MySource.GhorbaniCoeffs[3] * qPow(r_work, 2.0) / (1 + MySource.GhorbaniCoeffs[4] * qPow(r_work, 2) + MySource.GhorbaniCoeffs[5] * qPow(r_work, 4)) * 180 / 3.14;
    }


    PACurve->Working_point_dB_abs.x[0] = 20 * std::log10(r_work);
    PACurve->Working_point_dB_abs.y[0] = 10 * std::log10(qPow(A_work, 2)) + MySource.linear_gain_dB;

    PACurve->Working_point_linear_abs.x[0] = qPow(10.0, PACurve->Working_point_dB_abs.x[0] / 10.0);
    PACurve->Working_point_linear_abs.y[0] = qPow(10.0, PACurve->Working_point_dB_abs.y[0] / 10.0);

    PACurve->Working_point_dB_norm.x[0] = 20 * std::log10(r_work / r_sat);
    PACurve->Working_point_dB_norm.y[0] = 10 * std::log10(qPow(A_work, 2) / qPow(Asat, 2));

    PACurve->Working_point_linear_norm.x[0] = qPow(10.0, PACurve->Working_point_dB_norm.x[0] / 10.0);
    PACurve->Working_point_linear_norm.y[0] = qPow(10.0, PACurve->Working_point_dB_norm.y[0] / 10.0);

    if(FIR_enable) {
        ApplyFIRWithMemory(PACurve->P_out_abs.dB, PACurve->Phi, MySource.W_FIRCoeffs, 5);
        ApplyFIRWithMemory(PACurve->P_out_abs.linear, PACurve->Phi, MySource.W_FIRCoeffs, 5);
        ApplyFIRWithMemory(PACurve->P_out_norm.dB, PACurve->Phi, MySource.W_FIRCoeffs, 5);
        ApplyFIRWithMemory(PACurve->P_out_norm.linear, PACurve->Phi, MySource.W_FIRCoeffs, 5);
    }   
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
    MySource.SC_f_carrier = UISource.SC_f_carrier;
    MySource.SC_symrate = UISource.SC_symrate;
    MySource.SC_rolloff = UISource.SC_rolloff;
    MySource.SC_filter_length = UISource.SC_filter_length;
    MySource.SC_FilterType = UISource.SC_FilterType;

    MySource.OFDM_f_carrier = UISource.OFDM_f_carrier;
    MySource.OFDM_Nfft = UISource.OFDM_Nfft;
    MySource.OFDM_f_carrier = UISource.OFDM_f_carrier;
    MySource.OFDM_Nfft = UISource.OFDM_Nfft;
    MySource.OFDM_GB_DC = UISource.OFDM_GB_DC;
    MySource.OFDM_GB_Nyq = UISource.OFDM_GB_Nyq;
    MySource.OFDM_cycle_prefix = UISource.OFDM_cycle_prefix;

    MySource.FDMA_f_carrier = UISource.FDMA_f_carrier;
    MySource.FDMA_symrate = UISource.FDMA_symrate;
    MySource.FDMA_num_subcarriers = UISource.FDMA_num_subcarriers;
    MySource.FDMA_step_carrier = UISource.FDMA_step_carrier;

    MySource.oversampling = UISource.oversampling;
    MySource.fs = UISource.fs;
    MySource.SNRSig = UISource.SNRSig;
    MySource.bb_delta = UISource.bb_delta;

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

