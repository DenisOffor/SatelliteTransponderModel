#ifndef GRAPHPLOTTING_H
#define GRAPHPLOTTING_H

#include <QObject>
#include <QVector>
#include "qcustomplot.h"
#include <QScreen>
#include "HelpfullStructs.h"
#include "ofdm.h"

class GraphPlotting : public QObject
{
    Q_OBJECT
private:
    QVector<double> x, y;
    QVector<double> wp_x, wp_y;
    QCustomPlot* plotIdealSymConst;
    QCustomPlot* plotPaCurve;
    QCustomPlot* plotsOfConstellations[6];
    QCustomPlot* plotsOfTimeDomain[6];
    QCustomPlot* plotsOfPSD[6];
    void InitializePaCurvePlot(QWidget* GraphWidget);
    void InitializeIdealSymConstPlot(QWidget* GraphWidget);
    void InitializeTimeDomainPlotting(QVector<QWidget*> TimeDomainGraphWidgets);
    void InitializeConstellationsPlotting(QVector<QWidget*> ConstellationsGraphWidgets);
    void InitializePSDPlotting(QVector<QWidget*> PSDGraphWidgets);

public:
    //need to get GroupBox to constuctor in order to place QCustomPlot in widget in GropuBox
    GraphPlotting();
    void init(QVector<QWidget*> SetupGraphWidgets, QVector<QWidget*> ConstellationsGraphWidgets,
              QVector<QWidget*> TimeDomainGraphWidgets, QVector<QWidget*> PSDGraphWidgets);
    ~GraphPlotting();
    void PlotPaCurve(PaCurve& PaCurveData, const QList<QAction*> actions);
    void PlotConstellationsPlots(const Symbols& MySymbols);
    void PlotTimeDomainPlots(const OfdmResult& CurrentOfdm);
    void PlotPSDPlots(const QVector<QVector<double>>& PSDs, const QVector<QVector<double>>& freqs);

public slots:
    void PlotIdealSymConstellation(const QString ModType);
};

#endif // GRAPHPLOTTING_H
