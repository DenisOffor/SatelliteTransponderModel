#ifndef RESULTS_CSV_H
#define RESULTS_CSV_H
#pragma once

#include <QFile>
#include <QFileInfo>
#include <QDir>
#include <QTextStream>
#include <QString>
#include "HelpfullStructs.h"

inline QString csvNum(double x)
{
    return QString::number(x, 'g', 16);   // нормальная точность для MATLAB
}

inline QString resultsCsvHeader()
{
    return "IBO,"
           "OBO_noDPD,OBO_withDPD,"
           "BER_noDPD,BER_withDPD,"
           "EVM_noDPD,EVM_withDPD,"
           "PARP_noDPD,PARP_withDPD,"
           "P_formed_noDPD,P_formed_withDPD,"
           "P_emitted_noDPD,P_emitted_withDPD,"
           "Gain_noDPD,Gain_withDPD,"
           "ACLR_l_noDPD,ACLR_u_noDPD,"
           "ACLR_l_withDPD,ACLR_u_withDPD,"
           "NMSE_noDPD,NMSE_withDPD,"
           "BB\n";
}

// Вызывать 1 раз в начале новой серии измерений.
// Очищает старый файл и пишет заголовок.
inline bool initResultsCsv(const QString& filePath)
{
    QFileInfo fi(filePath);
    QDir().mkpath(fi.absolutePath());

    QFile file(filePath);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text | QIODevice::Truncate))
        return false;

    QTextStream out(&file);
    out << resultsCsvHeader();
    return true;
}

// Вызывать в конце КАЖДОГО цикла измерения.
// ibo передаешь отдельно, остальные поля берутся из ts.
template <typename TimeSignalT>
bool appendResultsCsv(const QString& filePath, double ibo, const TimeSignalT& ts)
{
    QFileInfo fi(filePath);
    QDir().mkpath(fi.absolutePath());

    bool needHeader = !fi.exists() || fi.size() == 0;

    QFile file(filePath);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text | QIODevice::Append))
        return false;

    QTextStream out(&file);
    out.setRealNumberNotation(QTextStream::SmartNotation);
    out.setRealNumberPrecision(16);

    if (needHeader)
        out << resultsCsvHeader();

    out << csvNum(ibo) << ","
        << csvNum(ts.CurrentOBO_noDPD) << "," << csvNum(ts.CurrentOBO_withDPD) << ","
        << csvNum(ts.BER_noDPD) << "," << csvNum(ts.BER_withDPD) << ","
        << csvNum(ts.EVM_noDPD) << "," << csvNum(ts.EVM_withDPD) << ","
        << csvNum(ts.PARP_noDPD) << "," << csvNum(ts.PARP_withDPD) << ","
        << csvNum(ts.P_formed_noDPD) << "," << csvNum(ts.P_formed_withDPD) << ","
        << csvNum(ts.P_emitted_noDPD) << "," << csvNum(ts.P_emitted_withDPD) << ","
        << csvNum(ts.Gain_noDPD) << "," << csvNum(ts.Gain_withDPD) << ","
        << csvNum(ts.ACPR_noDPD.first) << "," << csvNum(ts.ACPR_noDPD.second) << ","
        << csvNum(ts.ACPR_withDPD.first) << "," << csvNum(ts.ACPR_withDPD.second) << ","
        << csvNum(ts.NMSE_noDPD) << "," << csvNum(ts.NMSE_withDPD) << ","
        << csvNum(ts.BB)
        << "\n";

    return true;
}

#endif // RESULTS_CSV_H
