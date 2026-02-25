#ifndef GRAPHPLOTTING_H
#define GRAPHPLOTTING_H
#include <QVector>
#include <QObject>
#include "qcustomplot.h"
#include <QScreen>
#include "HelpfullStructs.h"
#include "ofdm.h"

class GraphPlotting : public QObject
{
    Q_OBJECT
private:
    QCustomPlot* plotIdealSymConst;
    QCustomPlot* plotPaCurve;
    QCustomPlot* plotsOfConstellations[6];
    QCustomPlot* plotsOfTimeDomain[6];
    QCustomPlot* plotsOfPSD[6];
    void InitializePaCurvePlot(QWidget* GraphWidget);
    void InitializeIdealSymConstPlot(QWidget* GraphWidget);
    void InitializeTimeDomainPlotting(std::vector<QWidget*> TimeDomainGraphWidgets);
    void InitializeConstellationsPlotting(std::vector<QWidget*> ConstellationsGraphWidgets);
    void InitializePSDPlotting(std::vector<QWidget*> PSDGraphWidgets);

public:
    //need to get GroupBox to constuctor in order to place QCustomPlot in widget in GropuBox
    GraphPlotting();
    void init(std::vector<QWidget*> SetupGraphWidgets, std::vector<QWidget*> ConstellationsGraphWidgets,
              std::vector<QWidget*> TimeDomainGraphWidgets, std::vector<QWidget*> PSDGraphWidgets);
    ~GraphPlotting();
    void PlotPaCurve(PaCurve& PaCurveData, const QList<QAction*> actions);
    void PlotConstellationsPlots(const Symbols& MySymbols);
    void PlotTimeDomainPlots(const GlobalResults& CurrentOfdm, bool rescale);
    void PlotPSDPlots(const std::vector<std::vector<double>>& PSDs, const std::vector<std::vector<double>>& freqs);

public slots:
    void PlotIdealSymConstellation(const QString ModType);
};

#endif // GRAPHPLOTTING_H
