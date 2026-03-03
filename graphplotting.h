#ifndef GRAPHPLOTTING_H
#define GRAPHPLOTTING_H
#include <QVector>
#include <QObject>
#include "qcustomplot.h"
#include <QScreen>
#include "HelpfullStructs.h"
#include "ofdm.h"
#include "./ui_mainwindow.h"

#define PA_CURVE_BIAS 0

class GraphPlotting : public QObject
{
    Q_OBJECT
private:
    Ui::MainWindow* Local_ui_copy;
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
    void changeGraphType();
    void saveGraphToFile();
    void copyGraphToClipboard();
    void onPlotDoubleClick(QMouseEvent *event);

public:
    //need to get GroupBox to constuctor in order to place QCustomPlot in widget in GropuBox
    GraphPlotting(Ui::MainWindow *ui, QObject *parent);
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
