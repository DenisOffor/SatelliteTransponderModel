#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QTreeWidgetItem>
#include <QListWidgetItem>
#include <QScreen>
#include <QDebug>
#include <QThread>
#include <QTimer>
#include "worker.h"
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
    void DPDTypeComboBoxTextChange(const QString& comboxString);
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
    QThread* thread;
    Worker* worker;
    QTimer *timer;

    void setupAdaptiveWindow();
    void setupSplitter();
    void setupPaCurveToolButton();
    void setupPACurvePlotting();
    void setupLabels();
    void setupMainPipelineTree();
    void SetupSelectedGraphsListWidget();
    void SetupSelectedDPDType();
    void SetupMainLogicWork();
    void SetupWorker();
    void SetupCyclePushBtn();
    void LockParChange();
    void UnLockParChange();

private slots:
    void handleResult();
    void cycleBtnClicked();
    void DPDRecalcBtnClicked();
    void CycleModeSlot();
signals:
    void startSimulation(SignalProcessing* sig, NeedToRecalc recalc);
};
#endif // MAINWINDOW_H
