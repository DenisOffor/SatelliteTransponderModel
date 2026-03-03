#include "signalprocessing.h"

SignalProcessing::SignalProcessing() : myFdma(mySC) {
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
    MySource = new Source;
    freq.resize(6);
    PSDs.resize(6);
}


SignalProcessing::~SignalProcessing() {
    delete PACurve;
    delete MySource;
}

void SignalProcessing::MainLogicWork(NeedToRecalc CurrentRecalcNeeds)
{
    QElapsedTimer timer;
    timer.start();
    if(CurrentRecalcNeeds.RecalcSymbols || CurrentRecalcNeeds.RecalcNoiseSym || CurrentRecalcNeeds.FullRecalc)
        GeneratePacksOfSymbols(CurrentRecalcNeeds);
    qDebug() << "Sym time:" << timer.elapsed() << "ms";

    timer.restart();
    if(CurrentRecalcNeeds.RecalcSig || CurrentRecalcNeeds.FullRecalc || CurrentRecalcNeeds.RecalcNoiseSig)
        TransmitSignalProcessing(CurrentRecalcNeeds);
    qDebug() << "Sig time:" << timer.elapsed() << "ms";

    timer.restart();
    if(CurrentRecalcNeeds.PARecalc || CurrentRecalcNeeds.FullRecalc)
        if(!CurrentRes.pa_sig.empty())
            PAProcessing(CurrentRecalcNeeds);
    qDebug() << "PA time:" << timer.elapsed() << "ms";

    timer.restart();
    if(!CurrentRes.pa_sig.empty())
        ReceiveSignalProcessing(CurrentRecalcNeeds);
    qDebug() << "Rec time:" << timer.elapsed() << "ms";

    timer.restart();
    MyMetricsEval.comparePSD(CurrentRes.tx_sig, CurrentRes.pa_sig, MySource->fs, MySource->oversampling, freq[0], PSDs[0], PSDs[1]);
    qDebug() << "PSD time:" << timer.elapsed() << "ms";
}

void SignalProcessing::GeneratePacksOfSymbols(NeedToRecalc& CurrentRecalcNeeds)
{
    if(CurrentRecalcNeeds.RecalcSymbols || CurrentRecalcNeeds.FullRecalc) {
        MySymbols.clear();
        if(MySource->SigType == "FDMA")
            for(int i = 0; i < MySource->FDMA_num_subcarriers; i++) {
                MySymbols.push_back(GenerateNSymbols());
                SymsAddNoise(MySymbols[i].tr_sym_clean, MySymbols[i].tr_sym_noisy);
            }
        else { MySymbols.push_back(GenerateNSymbols()); SymsAddNoise(MySymbols[0].tr_sym_clean, MySymbols[0].tr_sym_noisy); }
    }
    else {
        if(MySource->SigType == "FDMA")
            for(int i = 0; i < MySource->FDMA_num_subcarriers; i++)
                SymsAddNoise(MySymbols[i].tr_sym_clean, MySymbols[i].tr_sym_noisy);
        else { SymsAddNoise(MySymbols[0].tr_sym_clean, MySymbols[0].tr_sym_noisy); }
    }
    CurrentRecalcNeeds.RecalcSig = true;
}

Symbols SignalProcessing::GenerateNSymbols()
{
    Symbols temp;
    temp.tr_sym_clean.resize(MySource->NumSym);
    temp.tr_sym_noisy.resize(MySource->NumSym);
    temp.data.resize(MySource->NumSym);
    for (int i = 0; i < MySource->NumSym; i++)
        temp.data[i] = QRandomGenerator::global()->bounded(MySource->M);

    if(MySource->ModType == "BPSK") {
        for (int i = 0; i < MySource->NumSym; i++) {
            temp.tr_sym_clean[i] = std::polar(1.0, 2 * M_PI * temp.data[i] / MySource->M);
        }
    }
    else if(MySource->ModType == "QPSK" || MySource->ModType == "16QAM"
               || MySource->ModType == "64QAM") {
        int sqrtM = std::sqrt(MySource->M);
        double norm = std::sqrt(2.0 / 3.0 * (MySource->M - 1)); // UnitAveragePower

        for (int i = 0; i < MySource->NumSym; i++) {
            int s = temp.data[i];
            int I = 2 * (s % sqrtM) - sqrtM + 1;
            int Q = 2 * (s / sqrtM) - sqrtM + 1;
            temp.tr_sym_clean[i] = std::complex<double>(I, Q) / norm;
        }
    }
    else if(MySource->ModType == "APSK") {
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
    P_signal /= MySource->NumSym;

    double P_noise = P_signal / pow(10.0, MySource->SNRSymdB / 10.0);
    double sigma = std::sqrt(P_noise / 2);

    std::mt19937 rng(12345);
    std::normal_distribution<double> gauss(0.0, sigma);

    for (int i = 0; i < MySource->NumSym; i++) {
        temp.tr_sym_noisy[i] = temp.tr_sym_clean[i] +
                                     std::complex<double>(gauss(rng), gauss(rng));
    }

    return temp;
}

void SignalProcessing::SymsAddNoise(std::vector<std::complex<double>>& symbols_clean, std::vector<std::complex<double>>& symbols_noisy) {
    double P_signal = 0;
    for (auto& s : symbols_clean)
        P_signal += std::norm(s);
    P_signal /= MySource->NumSym;

    double P_noise = P_signal / pow(10.0, MySource->SNRSymdB / 10.0);
    double sigma = std::sqrt(P_noise / 2);

    std::mt19937 rng(12345);
    std::normal_distribution<double> gauss(0.0, sigma);

    for (int i = 0; i < MySource->NumSym; i++) {
        symbols_noisy[i] = symbols_clean[i] +
                            std::complex<double>(gauss(rng), gauss(rng));
    }
}


OfdmParams SignalProcessing::GetOfdmParams()
{
    OfdmParams ofdm_parms;
    ofdm_parms.CP = MySource->OFDM_cycle_prefix;
    ofdm_parms.GB_DC = MySource->OFDM_GB_DC;
    ofdm_parms.GB_Nyq = MySource->OFDM_GB_Nyq;
    ofdm_parms.Nfft = MySource->OFDM_Nfft;
    ofdm_parms.SNR_dB = MySource->SNRSig;
    ofdm_parms.fc = MySource->OFDM_f_carrier;
    ofdm_parms.fs = MySource->fs;
    ofdm_parms.oversampling = MySource->oversampling;
    return ofdm_parms;
}

ScParams SignalProcessing::GetSCParams()
{
    ScParams sc_params;
    sc_params.SC_FilterType = MySource->SC_FilterType;
    sc_params.SC_filter_length = MySource->SC_filter_length;
    sc_params.SC_f_carrier = MySource->SC_f_carrier;
    sc_params.SC_rolloff = MySource->SC_rolloff;
    sc_params.SC_symrate = MySource->SC_symrate;
    sc_params.SNR_dB = MySource->SNRSig;
    sc_params.fs = MySource->fs;
    sc_params.oversampling = MySource->oversampling;
    return sc_params;
}

FdmaParams SignalProcessing::GetFDMAParams()
{
    FdmaParams fdma_params;
    fdma_params.FDMA_f_carrier = MySource->FDMA_f_carrier;
    fdma_params.FDMA_symrate = MySource->FDMA_symrate;
    fdma_params.FDMA_num_subcarriers = MySource->FDMA_num_subcarriers;
    fdma_params.FDMA_step_carrier = MySource->FDMA_step_carrier;
    fdma_params.SNRSig = MySource->SNRSig;
    fdma_params.fs = MySource->fs;
    fdma_params.oversampling = MySource->oversampling;

    return fdma_params;
}


void SignalProcessing::TransmitSignalProcessing(NeedToRecalc& CurrentRecalcNeeds)
{
    if(MySource->SigType == "OFDM") {
        OfdmParams ofdm_parms = GetOfdmParams();
        if(CurrentRecalcNeeds.RecalcNoiseSig)
            myOfdm.changeAwgn(CurrentOfdmResults, ofdm_parms);
        else {
            CurrentOfdmResults = myOfdm.makeOfdm(MySymbols[0].tr_sym_noisy, ofdm_parms);
        }
        CurrentRes.clear();
        CurrentRes.resize(CurrentOfdmResults.tx.size());
        CurrentRes.tx_sig = CurrentOfdmResults.tx;
        CurrentRes.time = CurrentOfdmResults.t;
    }
    else if(MySource->SigType == "FDMA") {
        FdmaParams fdma_parms = GetFDMAParams();
        ScParams sc_params = GetSCParams();
        if(CurrentRecalcNeeds.RecalcNoiseSig)
            myFdma.changeAwgn(CurrentFdmaResults, fdma_parms);
        else
            CurrentFdmaResults = myFdma.generate(MySymbols, fdma_parms, sc_params);
        CurrentRes.clear();
        CurrentRes.resize(CurrentFdmaResults.tx.size());
        CurrentRes.tx_sig = CurrentFdmaResults.tx;
        CurrentRes.time = CurrentFdmaResults.t;
    }
    else if(MySource->SigType == "SC") {
        ScParams sc_params = GetSCParams();
        if(CurrentRecalcNeeds.RecalcNoiseSig && !CurrentRecalcNeeds.FullRecalc)
            mySC.changeAwgn(CurrentSCResults, sc_params);
        else {
            CurrentSCResults = mySC.makeSc(MySymbols[0].tr_sym_noisy, sc_params);
        }
        CurrentRes.clear();
        CurrentRes.resize(CurrentSCResults.tx.size());
        CurrentRes.tx_sig = CurrentSCResults.tx;
        CurrentRes.time = CurrentSCResults.t;
    }
    CurrentRecalcNeeds.PARecalc = true;
}

void SignalProcessing::PAProcessing(NeedToRecalc& CurrentRecalcNeeds)
{
    if(!CurrentRes.pa_sig.empty()) {
        CurrentRes.pa_sig = CurrentRes.tx_sig;
        if(MySource->PAModel == "Saleh")
            MyPAModels.SalehModel(CurrentRes.pa_sig, MySource->SalehCoeffs, MySource->linear_gain_dB, MySource->IBO_dB);
        else if (MySource->PAModel == "Rapp")
            MyPAModels.RappModel(CurrentRes.pa_sig, MySource->RappCoeffs, MySource->linear_gain_dB, MySource->IBO_dB);
        else if (MySource->PAModel == "Ghorbani")
            MyPAModels.GhorbaniModel(CurrentRes.pa_sig, MySource->GhorbaniCoeffs, MySource->linear_gain_dB, MySource->IBO_dB);
        else if (MySource->PAModel == "Wiener") {
            if(MySource->StaticNonlinModel == "Saleh")
                MyPAModels.WienerModel(CurrentRes.pa_sig, MySource->StaticNonlinModel,
                                       MySource->SalehCoeffs, MySource->FIRCoeffs, MySource->linear_gain_dB, MySource->IBO_dB);
            else if(MySource->StaticNonlinModel == "Rapp")
                MyPAModels.WienerModel(CurrentRes.pa_sig, MySource->StaticNonlinModel,
                                       MySource->RappCoeffs, MySource->FIRCoeffs, MySource->linear_gain_dB, MySource->IBO_dB);
            else if(MySource->StaticNonlinModel == "Ghorbani")
                MyPAModels.WienerModel(CurrentRes.pa_sig, MySource->StaticNonlinModel,
                                       MySource->GhorbaniCoeffs, MySource->FIRCoeffs, MySource->linear_gain_dB, MySource->IBO_dB);
        }
    }
}

void SignalProcessing::ReceiveSignalProcessing(NeedToRecalc CurrentRecalcNeeds)
{
    OfdmParams ofdm_parms = GetOfdmParams();
    ScParams sc_params = GetSCParams();
    FdmaParams fdma_params = GetFDMAParams();
    if(MySource->SigType == "OFDM")
        MySymbols[0].rec_sym_noisy = myOfdm.ofdm_demodulate(CurrentRes.pa_sig, MySymbols[0].tr_sym_clean, ofdm_parms);
    else if(MySource->SigType == "SC")
        MySymbols[0].rec_sym_noisy = mySC.demodulateSignal(CurrentRes.pa_sig, sc_params, sc_params);
    else if(MySource->SigType == "FDMA") {
        std::vector<std::vector<std::complex<double>>> temp = myFdma.demodulate(CurrentRes.pa_sig, fdma_params, sc_params);
        for(int i = 0; i < MySymbols.size(); ++i)
            MySymbols[i].rec_sym_noisy = temp[i];
    }
}

void SignalProcessing::CalcPaCurve(const QString model_type, const std::vector<double> Coefs, const int IBO_dB, const int LinearGaindB, bool FIR_enable,
                                        const std::vector<double> FIR_Coefs)
{
    double linear_gain = qPow(10.0, LinearGaindB / 20);
    if(model_type == "Saleh") {
        for(int i = 0; i < PACurve->point_num; i ++) {
            // AM/AM (unnormalized)
            PACurve->Aout[i] = (Coefs[0] * PACurve->r_in[i]) / (1 + Coefs[1] * qPow(PACurve->r_in[i], 2.0)) ;
            // AM/PM
            PACurve->Phi[i] = (Coefs[2] * qPow(PACurve->r_in[i], 2)) /  (1 + Coefs[3] * qPow(PACurve->r_in[i], 2)) * 180 / 3.14;
        }
    }
    else if (model_type == "Rapp") {
        for(int i = 0; i < PACurve->point_num; i ++) {
            // AM/AM (unnormalized)
            PACurve->Aout[i] = linear_gain * PACurve->r_in[i] / qPow(1 + qPow(PACurve->r_in[i] / Coefs[0], 2 * Coefs[1]), 1 / (2 * Coefs[1]));
            // AM/PM
            PACurve->Phi[i] = 0;
        }
    }
    else if (model_type == "Ghorbani") {
        for(int i = 0; i < PACurve->point_num; i ++) {
            // AM/AM (unnormalized)
            PACurve->Aout[i] = Coefs[0] * PACurve->r_in[i] / (1 + Coefs[1] * qPow(PACurve->r_in[i], 2) + Coefs[2] * qPow(PACurve->r_in[i], 4));
            // AM/PM
            PACurve->Phi[i] = Coefs[3] * qPow(PACurve->r_in[i], 2) / (1 + Coefs[4] * qPow(PACurve->r_in[i], 2) + Coefs[5] * qPow(PACurve->r_in[i], 4)) * 180 / 3.14;
        }
    }

    // Output power in dB and linear
    for(int i = 0; i < PACurve->point_num; i ++) {
        // P out in dB with linear gain
        PACurve->P_out_abs.dB[i] =  10 * std::log10(qPow(PACurve->Aout[i], 2)) + LinearGaindB;
        PACurve->P_out_abs.linear[i] = qPow(10.0, PACurve->P_out_abs.dB[i] / 10.0);
    }

    // Make artificial saturation
    // --- ЭМПИРИЧЕСКАЯ сатурация (по факту, без аналитической асимптоты) ---
    auto maxIt = std::max_element(PACurve->P_out_abs.dB.begin(), PACurve->P_out_abs.dB.end());
    double Psat = *maxIt;               // значение максимума
    int maxIndex = maxIt - PACurve->P_out_abs.dB.begin();    // индекс максимума
    double P_in_sat_dB = PACurve->P_in_abs.dB[maxIndex];
    double P_in_sat_linear = PACurve->P_in_abs.linear[maxIndex];
    double r_sat = PACurve->r_in[maxIndex];
    // мощность насыщения (включая gain_linear) — используется для нормировки дБ
    double P_sat_dB_with_gain = Psat + LinearGaindB;

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
    double IBO_linear = qPow(10.0, IBO_dB / 10.0);
    double r_work = r_sat / std::sqrt(IBO_linear);

    double A_work;
    if(model_type == "Saleh") {
        A_work = (Coefs[0] * r_work) / (1 + Coefs[1] * qPow(r_work, 2));
        PACurve->Phi_work_grad[0] = (Coefs[2] * qPow(r_work, 2)) / (1 + Coefs[3] * qPow(r_work, 2)) * 180 / 3.14;
    }
    else if(model_type == "Rapp") {
        A_work = r_work / qPow(1 + qPow(r_work / Coefs[0], 2 * Coefs[1]), 1 / (2 * Coefs[1]));
        PACurve->Phi_work_grad[0] = 0;
    }
    else if(model_type == "Ghorbani") {
        A_work = Coefs[0] * r_work / (1 + Coefs[1] * qPow(r_work, 2) + Coefs[2] * qPow(r_work, 4));
        PACurve->Phi_work_grad[0] = Coefs[3] * qPow(r_work, 2.0) / (1 + Coefs[4] * qPow(r_work, 2) + Coefs[5] * qPow(r_work, 4)) * 180 / 3.14;
    }


    PACurve->Working_point_dB_abs.x[0] = 20 * std::log10(r_work);
    PACurve->Working_point_dB_abs.y[0] = 10 * std::log10(qPow(A_work, 2)) + LinearGaindB;

    PACurve->Working_point_linear_abs.x[0] = qPow(10.0, PACurve->Working_point_dB_abs.x[0] / 10.0);
    PACurve->Working_point_linear_abs.y[0] = qPow(10.0, PACurve->Working_point_dB_abs.y[0] / 10.0);

    PACurve->Working_point_dB_norm.x[0] = 20 * std::log10(r_work / r_sat);
    PACurve->Working_point_dB_norm.y[0] = 10 * std::log10(qPow(A_work, 2)) - P_sat_dB_with_gain;

    PACurve->Working_point_linear_norm.x[0] = qPow(10.0, PACurve->Working_point_dB_norm.x[0] / 10.0);
    PACurve->Working_point_linear_norm.y[0] = qPow(10.0, PACurve->Working_point_dB_norm.y[0] / 10.0);

    if(FIR_enable) {
        ApplyFIRWithMemory(PACurve->P_out_abs.dB, PACurve->Phi, FIR_Coefs, 3);
        ApplyFIRWithMemory(PACurve->P_out_abs.linear, PACurve->Phi, FIR_Coefs, 3);
        ApplyFIRWithMemory(PACurve->P_out_norm.dB, PACurve->Phi, FIR_Coefs, 3);
        ApplyFIRWithMemory(PACurve->P_out_norm.linear, PACurve->Phi, FIR_Coefs, 3);
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
        phase[n] = std::arg(output[n]);
    }
}

void SignalProcessing::DataUpdate(Source &UISource)
{
    MySource->ModType = UISource.ModType;
    MySource->SNRSymdB = UISource.SNRSymdB;
    MySource->NumSym = UISource.NumSym;
    MySource->M = UISource.M;
    MySource->SigType = UISource.SigType;

    MySource->SC_f_carrier = UISource.SC_f_carrier;
    MySource->SC_symrate = UISource.SC_symrate;
    MySource->SC_rolloff = UISource.SC_rolloff;
    MySource->SC_filter_length = UISource.SC_filter_length;
    MySource->SC_FilterType = UISource.SC_FilterType;

    MySource->OFDM_f_carrier = UISource.OFDM_f_carrier;
    MySource->OFDM_Nfft = UISource.OFDM_Nfft;
    MySource->OFDM_f_carrier = UISource.OFDM_f_carrier;
    MySource->OFDM_Nfft = UISource.OFDM_Nfft;
    MySource->OFDM_GB_DC = UISource.OFDM_GB_DC;
    MySource->OFDM_GB_Nyq = UISource.OFDM_GB_Nyq;
    MySource->OFDM_cycle_prefix = UISource.OFDM_cycle_prefix;

    MySource->FDMA_f_carrier = UISource.FDMA_f_carrier;
    MySource->FDMA_symrate = UISource.FDMA_symrate;
    MySource->FDMA_num_subcarriers = UISource.FDMA_num_subcarriers;
    MySource->FDMA_step_carrier = UISource.FDMA_step_carrier;

    MySource->oversampling = UISource.oversampling;
    MySource->fs = UISource.fs;
    MySource->SNRSig = UISource.SNRSig;

    MySource->PAModel = UISource.PAModel;
    MySource->IBO_dB = UISource.IBO_dB;
    MySource->linear_gain_dB = UISource.linear_gain_dB;
    MySource->SalehCoeffs = UISource.SalehCoeffs;
    MySource->RappCoeffs = UISource.RappCoeffs;
    MySource->GhorbaniCoeffs = UISource.GhorbaniCoeffs;
    MySource->FIRCoeffs = UISource.FIRCoeffs;
    MySource->StaticNonlinModel = UISource.StaticNonlinModel;
}

Symbols &SignalProcessing::getSymbols()
{
        return MySymbols[0];
}

GlobalResults &SignalProcessing::getTimeSignal()
{
        return CurrentRes;
}

std::vector<std::vector<double> > SignalProcessing::getFreq()
{
    return freq;
}

std::vector<std::vector<double> > SignalProcessing::getPSDs()
{
    return PSDs;
}

PaCurve &SignalProcessing::getPaCurve()
{
    return *PACurve;
}








