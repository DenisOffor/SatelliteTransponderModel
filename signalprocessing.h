#ifndef SIGNALPROCESSING_H
#define SIGNALPROCESSING_H
#include <QVector>
#include "HelpfullStructs.h"
#include "ofdm.h"
#include "pamodels.h"
#include "metricseval.h"

class SignalProcessing
{
private:
    PaCurve* PACurve;
    QVector<Symbols> MySymbols;
    Source* MySource;
    PAModels MyPAModels;
    OFDM myOfdm;
    MetricsEval MyMetricsEval;
    OfdmResult CurrentOfdmResults;
    GlobalResults CurrentRes;
    QVector<QVector<double>> freq;
    QVector<QVector<double>> PSDs;

    //functions
    Symbols GenerateNSymbols();
    OfdmParams GetOfdmParams();
    void TransmitSignalProcessing(NeedToRecalc CurrentRecalcNeeds);
    void PAProcessing(NeedToRecalc CurrentRecalcNeeds);
    void ReceiveSignalProcessing(NeedToRecalc CurrentRecalcNeeds);

public:
    SignalProcessing();
    ~SignalProcessing();
    void GeneratePacksOfSymbols();
    PaCurve& CalcPaCurve(const QString model_type, const QVector<double> SalehCoefs, const int IBO_dB, const int LinearGain);
    void DataUpdate(Source& UISource);
    Symbols& getSymbols();
    OfdmResult& getTimeSignal();
    QVector<QVector<double>> getFreq();
    QVector<QVector<double>> getPSDs();

public slots:
    void MainLogicWork(NeedToRecalc CurrentRecalcNeeds);
};

#endif // SIGNALPROCESSING_H
