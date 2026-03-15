#include "worker.h"

void Worker::process(SignalProcessing* MySigProc, NeedToRecalc CurrentRecalcNeeds)
{
    MySigProc->MainLogicWork(CurrentRecalcNeeds);
    emit resultReady();
}
