#ifndef SIGNALPROCESSING_H
#define SIGNALPROCESSING_H
#include <QVector>
#include "HelpfullStructs.h"
#include "ofdm.h"

class SignalProcessing
{
private:
    PaCurve* PACurve;
    QVector<Symbols> MySymbols;
    Source* MySource;
    OFDM myOfdm;
    OfdmResult CurrentOfdmResults;

    //functions
    Symbols GenerateNSymbols();
    OfdmParams GetOfdmParams();

public:
    SignalProcessing();
    ~SignalProcessing();
    void GeneratePacksOfSymbols();
    PaCurve& CalcPaCurve(const QString model_type, const QVector<double> SalehCoefs, const int IBO_dB, const int LinearGain);
    void DataUpdate(Source& UISource);
    Symbols& getSymbols();
    OfdmResult& getTimeSignal();

public slots:
    void MainLogicWork(NeedToRecalc CurrentRecalcNeeds);
};

#endif // SIGNALPROCESSING_H
