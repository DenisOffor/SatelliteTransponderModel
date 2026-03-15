#ifndef WORKER_H
#define WORKER_H

#include <QObject>
#include <signalprocessing.h>

class Worker : public QObject
{
    Q_OBJECT
public slots:
    void process(SignalProcessing* sig, NeedToRecalc recalc);

signals:
    void resultReady();
};

#endif // WORKER_H
