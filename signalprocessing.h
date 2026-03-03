#ifndef SIGNALPROCESSING_H
#define SIGNALPROCESSING_H
#include "HelpfullStructs.h"
#include "ofdm.h"
#include "sc.h"
#include "fdma.h"
#include "pamodels.h"
#include "metricseval.h"
#include "QElapsedTimer"

class SignalProcessing
{
private:
    PaCurve* PACurve;
    std::vector<Symbols> MySymbols;
    Source* MySource;
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

    //functions
    Symbols GenerateNSymbols();
    OfdmParams GetOfdmParams();
    ScParams GetSCParams();
    FdmaParams GetFDMAParams();
    void TransmitSignalProcessing(NeedToRecalc& CurrentRecalcNeeds);
    void PAProcessing(NeedToRecalc& CurrentRecalcNeeds);
    void ReceiveSignalProcessing(NeedToRecalc CurrentRecalcNeeds);
    void SymsAddNoise(std::vector<std::complex<double>>& symbols_clean,
                      std::vector<std::complex<double>>& symbols_noisy);
public:
    SignalProcessing();
    ~SignalProcessing();
    void GeneratePacksOfSymbols(NeedToRecalc& CurrentRecalcNeeds);
    void CalcPaCurve(const QString model_type, const std::vector<double> SalehCoefs, const int IBO_dB, const int LinearGain, bool FIR_enable,
                        const std::vector<double> FIR_Coefs);
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
