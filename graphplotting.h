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


struct GraphActions {
    QList<QAction*> common;     // Общие для всех графиков
    QList<QAction*> specific;   // Специфичные для конкретного типа
};

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

    QMap<QString, GraphActions> m_graphActions;  // Тип графика -> действия
    QString m_currentGraphType;
    GraphActions CommonGraphActions;
    GraphActions IdealSymConstGraphActions;
    GraphActions PaCurveGraphActions;
    GraphActions SymConstGraphActions;
    GraphActions TimeDomainGraphActions;
    GraphActions PSDGraphActions;

    void InitializePaCurvePlot(QWidget* GraphWidget);
    void InitializeIdealSymConstPlot(QWidget* GraphWidget);
    void InitializeTimeDomainPlotting(std::vector<QWidget*> TimeDomainGraphWidgets);
    void InitializeConstellationsPlotting(std::vector<QWidget*> ConstellationsGraphWidgets);
    void InitializePSDPlotting(std::vector<QWidget*> PSDGraphWidgets);

    void setupGraphActions();
    void updateMenuForGraphType(QMenu &menu, const QString &graphType);
    void changeGraphType();
    void saveGraphToFile();
    void copyGraphToClipboard();
    void onPlotClick(QMouseEvent *event);

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
