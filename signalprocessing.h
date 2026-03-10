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
    Source* MySource;
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
    OfdmParams GetOfdmParams();
    ScParams GetSCParams();
    FdmaParams GetFDMAParams();
    void TransmitSignalProcessing(Source& source, NeedToRecalc& CurrentRecalcNeeds);
    void PAProcessing(Source& source, NeedToRecalc& CurrentRecalcNeeds);
    void ReceiveSignalProcessing(Source& source, NeedToRecalc CurrentRecalcNeeds);
    void SymsAddNoise(std::vector<std::complex<double>>& symbols_clean,
                      std::vector<std::complex<double>>& symbols_noisy);
    int DemodulateSymbol(const std::complex<double>& r, const std::vector<std::complex<double>>& constellation);
    void Demodulate(Symbols& symbols, const std::vector<std::complex<double>>& constellation);
    void InitializeConstellations();
    void RecalcDPD(NeedToRecalc& CurrentRecalcNeeds);
    std::vector<std::complex<double>>& getCurrentConstellation();
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
