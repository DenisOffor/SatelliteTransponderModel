#ifndef SIGNALPROCESSING_H
#define SIGNALPROCESSING_H
#include "HelpfullStructs.h"
#include "ofdm.h"
#include "sc.h"
#include "fdma.h"
#include "pamodels.h"
#include "metricseval.h"
#include "QElapsedTimer"
#include "dpd.h"

class SignalProcessing
{
private:
    PaCurve* PACurve;
    std::vector<Symbols> MySymbols;
    Source MySource;
    Source TrainSource;
    DPD mydpd;
    PAModels MyPAModels;
    OFDM myOfdm;
    SC mySC;
    FDMA myFdma;
    OfdmResult CurrentOfdmResults;
    ScResult CurrentSCResults;
    FdmaResult CurrentFdmaResults;
    GlobalResults CurrentRes;
    MetricsEval MyMetricsEval;
    std::vector<std::vector<double>> freq;
    std::vector<std::vector<double>> PSDs;

    std::vector<std::complex<double>> BPSK_const;
    std::vector<std::complex<double>> QPSK_const;
    std::vector<std::complex<double>> QAM16_const;
    std::vector<std::complex<double>> QAM64_const;

    //functions
    Symbols GenerateNSymbols(Source& source);
    OfdmParams GetOfdmParams(Source& source);
    ScParams GetSCParams(Source& source);
    FdmaParams GetFDMAParams(Source& source);
    void TransmitSignalProcessing(Source& source, std::vector<Symbols>& symbols, NeedToRecalc& CurrentRecalcNeeds, GlobalResults& CurRes);
    void PAProcessing(Source& source, NeedToRecalc& CurrentRecalcNeeds, GlobalResults& CurRes);
    void ReceiveSignalProcessing(Source& source, std::vector<Symbols>& symbols, NeedToRecalc CurrentRecalcNeeds, GlobalResults& CurRes);
    void SymsAddNoise(Source& source, std::vector<std::complex<double>>& symbols_clean,
                      std::vector<std::complex<double>>& symbols_noisy);
    int DemodulateSymbol(const std::complex<double>& r, const std::vector<std::complex<double>>& constellation);
    void Demodulate(Symbols& symbols, const std::vector<std::complex<double>>& constellation);
    void InitializeConstellations();
    void RecalcDPD(NeedToRecalc& CurrentRecalcNeeds);
    std::vector<std::complex<double>>& getCurrentConstellation(Source& source);
public:
    SignalProcessing();
    ~SignalProcessing();
    void GeneratePacksOfSymbols(std::vector<Symbols>& Symbols, Source& source, NeedToRecalc& CurrentRecalcNeeds);
    void CalcPaCurve();
    void ApplyFIRWithMemory(std::vector<double>& amplitude, std::vector<double>& phase, const std::vector<double>& FIR_Coefs, int numTaps);
    void DataUpdate(Source& UISource);
    Symbols& getSymbols();
    GlobalResults& getTimeSignal();
    std::vector<std::vector<double>> getFreq();
    std::vector<std::vector<double>> getPSDs();
    PaCurve& getPaCurve();

public slots:
    void MainLogicWork(NeedToRecalc CurrentRecalcNeeds);
};

#endif // SIGNALPROCESSING_H
