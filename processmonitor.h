#ifndef PROCESSMONITOR_H
#define PROCESSMONITOR_H

#pragma once

#include <QObject>
#include <QString>

class ProcessMonitor : public QObject
{
    Q_OBJECT
public:
    explicit ProcessMonitor(QObject *parent = nullptr);

    void reset();
    void update();

    double cpuUsagePercent() const;

    double workingSetMB() const;
    double privateUsageMB() const;
    double peakWorkingSetMB() const;

    QString cpuText(int precision = 1) const;
    QString ramText(int precision = 1) const;

private:
    quint64 m_prevProcTime100ns = 0;
    quint64 m_prevWallTime100ns = 0;
    bool m_hasPrevSample = false;

    double m_cpuUsagePercent = 0.0;

    double m_workingSetMB = 0.0;
    double m_privateUsageMB = 0.0;
    double m_peakWorkingSetMB = 0.0;
};
#endif // PROCESSMONITOR_H
