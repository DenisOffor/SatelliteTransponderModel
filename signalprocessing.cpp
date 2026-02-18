#include "signalprocessing.h"

SignalProcessing::SignalProcessing() {
    PACurve = new PaCurve(200);
    // Input power in dB and linear
    double dB_start = -20;
    double dB_end = 5;
    double dB_step = (dB_end - dB_start) / (PACurve->point_num - 1);
    for(int i = 0; i < PACurve->point_num; i ++) {
        PACurve->P_in_abs.dB[i] = dB_start + i * dB_step;
        PACurve->P_in_abs.linear[i] = qPow(10.0, PACurve->P_in_abs.dB[i] / 10.0);
        PACurve->r_in[i] = qSqrt(PACurve->P_in_abs.linear[i]);
    }
    MySource = new Source;
}


SignalProcessing::~SignalProcessing() {
    delete PACurve;
    delete MySource;
}

void SignalProcessing::MainLogicWork(NeedToRecalc CurrentRecalcNeeds)
{
    GeneratePacksOfSymbols();
    if(MySource->SigType == "OFDM") {        
        OfdmParams ofdm_parms = GetOfdmParams();
        CurrentOfdmResults = myOfdm.makeOfdm(MySymbols[0].tr_sym_noisy, ofdm_parms);
        MySymbols[0].rec_sym_noisy = myOfdm.ofdm_demodulate(CurrentOfdmResults.tx, MySymbols[0].tr_sym_clean, ofdm_parms);
    }
    else if(MySource->SigType == "FDMA") {

    }
    else if(MySource->SigType == "SC")
    return;
}

void SignalProcessing::GeneratePacksOfSymbols()
{
    MySymbols.clear();
    if(MySource->SigType == "FDMA")
        for(int i = 0; i < MySource->FDMA_num_subcarriers; i++)
            MySymbols.append(GenerateNSymbols());
    else MySymbols.append(GenerateNSymbols());
}

Symbols SignalProcessing::GenerateNSymbols()
{
    Symbols temp;
    temp.resize(MySource->NumSym);
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

PaCurve& SignalProcessing::CalcPaCurve(const QString model_type, const QVector<double> Coefs, const int IBO_dB, const int LinearGaindB)
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
    else if (model_type == "Memory Polynimial") {}

    // Output power in dB and linear
    for(int i = 0; i < PACurve->point_num; i ++) {
        // P out in dB with linear gain
        PACurve->P_out_abs.dB[i] =  10 * std::log10(qPow(PACurve->Aout[i], 2)) + LinearGaindB;
        PACurve->P_out_abs.linear[i] = qPow(10.0, PACurve->P_out_abs.dB[i] / 10.0);
    }

    // Make artificial saturation
    // --- ЭМПИРИЧЕСКАЯ сатурация (по факту, без аналитической асимптоты) ---
    auto maxIt = std::max_element(PACurve->P_out_abs.dB.begin(), PACurve->P_out_abs.dB.end());
    double Asat = *maxIt;               // значение максимума
    int maxIndex = maxIt - PACurve->P_out_abs.dB.begin();    // индекс максимума
    double P_in_sat_dB = PACurve->P_in_abs.dB[maxIndex];
    double P_in_sat_linear = PACurve->P_in_abs.linear[maxIndex];
    double r_sat = PACurve->r_in[maxIndex];
    // мощность насыщения (включая gain_linear) — используется для нормировки дБ
    double P_sat_dB_with_gain = 10 * std::log10(qPow(Asat, 2)) + LinearGaindB;

    //Нормированная выходная мощность: вычитаем Psat (gain_linear автоматически уйдёт)
    for (int i = 0; i < PACurve->point_num; i++) {
        PACurve->P_out_norm.dB[i] = PACurve->P_out_abs.dB[i] - P_sat_dB_with_gain;
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
        A_work = linear_gain * r_work / qPow(1 + qPow(r_work / Coefs[0], 2 * Coefs[1]), 1 / (2 * Coefs[1]));
        PACurve->Phi_work_grad[0] = 0;
    }
    else if(model_type == "Ghorbani") {
        A_work = Coefs[0] * r_work / (1 + Coefs[1] * qPow(r_work, 2) + Coefs[2] * qPow(r_work, 4));
        PACurve->Phi_work_grad[0] = Coefs[3] * qPow(r_work, 2.0) / (1 + Coefs[4] * qPow(r_work, 2) + Coefs[5] * qPow(r_work, 4)) * 180 / 3.14;
    }
    else if (model_type == "Memory Polynomial") {}


    PACurve->Working_point_dB_abs.x[0] = 20 * std::log10(r_work);
    PACurve->Working_point_dB_abs.y[0] = 10 * std::log10(qPow(A_work, 2)) + LinearGaindB;

    PACurve->Working_point_linear_abs.x[0] = qPow(10.0, PACurve->Working_point_dB_abs.x[0] / 10.0);
    PACurve->Working_point_linear_abs.y[0] = qPow(10.0, PACurve->Working_point_dB_abs.y[0] / 10.0);

    PACurve->Working_point_dB_norm.x[0] = 20 * std::log10(r_work / r_sat);
    PACurve->Working_point_dB_norm.y[0] = 10 * std::log10(qPow(A_work, 2) / qPow(Asat, 2));

    PACurve->Working_point_linear_norm.x[0] = r_work / r_sat;
    PACurve->Working_point_linear_norm.y[0] = qPow(A_work, 2) / qPow(Asat, 2);

    return *PACurve;
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
}

Symbols &SignalProcessing::getSymbols()
{
    if(!MySymbols.empty())
        return MySymbols[0];
}

OfdmResult &SignalProcessing::getTimeSignal()
{
        return CurrentOfdmResults;
}










