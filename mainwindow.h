#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QTreeWidgetItem>
#include <QListWidgetItem>
#include <QScreen>
#include <QDebug>
#include "graphplotting.h"
#include "signalprocessing.h"

QT_BEGIN_NAMESPACE
namespace Ui {
class MainWindow;
}
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
protected:
public slots:
    void onPipelineItemChanged(QTreeWidgetItem* current, QTreeWidgetItem*);
    void SignalTypeComboBoxTextChanged(const QString& comboxString);
    void PaModelTypeComboBoxTextChanged(const QString& comboxString);
    void onGraphsListItemChanged(QListWidgetItem* current, QListWidgetItem*);
    void PaCurvePlot();
    void MakeMainCalcAndPlot();
    void DataUpdate();
    void FirstDataUpdate();

private:
    Ui::MainWindow *ui;
    GraphPlotting Graphs;
    SignalProcessing MySigProc;
    Source UISource;
    NeedToRecalc CurrentRecalcNeeds;

    void setupAdaptiveWindow();
    void setupSplitter();
    void setupPaCurveToolButton();
    void setupPACurvePlotting();
    void setupLabels();
    void setupMainPipelineTree();
    void SetupSelectedGraphsListWidget();
    void SetupMainLogicWork();
};
#endif // MAINWINDOW_H
