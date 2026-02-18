#include "mainwindow.h"
#include "./ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow),
    Graphs(),
    MySigProc()
{
    ui->setupUi(this);
    QVector<QWidget*> SetupGraphWidgets;
    SetupGraphWidgets << ui->SymConstGroupBoxWidget << ui->PAPlotCurve_Widget;
    QVector<QWidget*> ConstellationsGraphWidgets;
    ConstellationsGraphWidgets << ui->ConstelPageGraph1_Widget << ui->ConstelPageGraph2_Widget << ui->ConstelPageGraph3_Widget
                 << ui->ConstelPageGraph4_Widget << ui->ConstelPageGraph5_Widget << ui->ConstelPageGraph6_Widget;
    QVector<QWidget*> TimeDomainGraphWidgets;
    TimeDomainGraphWidgets << ui->TimeDomainGraph1_Widget << ui->TimeDomainGraph2_Widget << ui->TimeDomainGraph3_Widget
                               << ui->TimeDomainGraph4_Widget << ui->TimeDomainGraph5_Widget << ui->TimeDomainGraph6_Widget;
    Graphs.init(SetupGraphWidgets, ConstellationsGraphWidgets, TimeDomainGraphWidgets);

    //Connecting chosen signal type with settings of this signal and setting starting page
    connect(ui->SignalTypeComboBox, &QComboBox::currentTextChanged, this, &MainWindow::SignalTypeComboBoxTextChanged);
    ui->PagesOfSigParametersStack->setCurrentWidget(ui->OFDM_page);

    //Connecting type of modulation and his signal constellation
    connect(ui->ModTypeComboBox, &QComboBox::currentTextChanged, &Graphs, &GraphPlotting::PlotIdealSymConstellation);

    FirstDataUpdate();
    setupMainPipelineTree();
    setupLabels();
    setupAdaptiveWindow();
    setupSplitter();
    setupPaCurveToolButton();
    setupPACurvePlotting();
    SetupSelectedGraphsListWidget();
    SetupMainLogicWork();
    CurrentRecalcNeeds.init();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::MakeMainCalcAndPlot()
{
    QObject *senderObj = sender();
    MySigProc.MainLogicWork(CurrentRecalcNeeds);
    Graphs.PlotConstellationsPlots(MySigProc.getSymbols());
    Graphs.PlotTimeDomainPlots(MySigProc.getTimeSignal());
}

void MainWindow::PaCurvePlot() {
    QVector<double> Coeffs;
    QString model_type = ui->PAModel_ComboBox->currentText();
    if(model_type == "Saleh") {
        Coeffs << ui->SalehCoef1_doubleSpinBox->value() << ui->SalehCoef2_doubleSpinBox->value()
               << ui->SalehCoef3_doubleSpinBox->value() << ui->SalehCoef4_doubleSpinBox->value();
        Graphs.PlotPaCurve(MySigProc.CalcPaCurve(model_type, Coeffs, ui->IBO_SpinBox->value(),
                                                 ui->LinearGain_SpinBox->value()),  ui->PACurveToolButton->menu()->actions());
    }
    else if (model_type == "Rapp") {
        Coeffs << ui->RappAsatCoef_doubleSpinBox->value() << ui->RappPCoef_doubleSpinBox->value();
        Graphs.PlotPaCurve(MySigProc.CalcPaCurve(model_type, Coeffs, ui->IBO_SpinBox->value(),
                                                   ui->LinearGain_SpinBox->value()), ui->PACurveToolButton->menu()->actions());
    }
    else if (model_type == "Ghorbani") {
        Coeffs << ui->GhorbaniCoef1_doubleSpinBox->value() << ui->GhorbaniCoef2_doubleSpinBox->value()
               << ui->GhorbaniCoef3_doubleSpinBox->value() << ui->GhorbaniCoef4_doubleSpinBox->value()
               << ui->GhorbaniCoef5_doubleSpinBox->value() << ui->GhorbaniCoef6_doubleSpinBox->value();
        Graphs.PlotPaCurve(MySigProc.CalcPaCurve(model_type, Coeffs, ui->IBO_SpinBox->value(),
                                                   ui->LinearGain_SpinBox->value()), ui->PACurveToolButton->menu()->actions());
    }
    else if (model_type == "Memory Polynomial") {}
}

void MainWindow::DataUpdate()
{
    if(ui->ModTypeComboBox->currentText() == "BPSK") UISource.M = 2;
    else if(ui->ModTypeComboBox->currentText() == "QPSK") UISource.M = 4;
    else if(ui->ModTypeComboBox->currentText() == "16QAM") UISource.M = 16;
    else if(ui->ModTypeComboBox->currentText() == "64QAM") UISource.M = 64;

    QObject *senderObj = sender();
    CurrentRecalcNeeds.clear();
    if(senderObj == ui->ModTypeComboBox) { UISource.ModType = ui->ModTypeComboBox->currentText(); CurrentRecalcNeeds.FullRecalc = true; }
    if(senderObj == ui->SNRSymSpinBox) { UISource.SNRSymdB = ui->SNRSymSpinBox->value(); CurrentRecalcNeeds.RecalcNoiseSym = true;
                                         CurrentRecalcNeeds.RecalcOfdmSig = true; CurrentRecalcNeeds.RecalcNoiseSig = true; }
    if(senderObj == ui->NumSymSpinBox) { UISource.NumSym = ui->NumSymSpinBox->value(); CurrentRecalcNeeds.FullRecalc = true; }
    if(senderObj == ui->SignalTypeComboBox) { UISource.SigType = ui->SignalTypeComboBox->currentText(); CurrentRecalcNeeds.FullRecalc = true; }

    if(senderObj == ui->SC_fc_SpinBox) { UISource.SC_f_carrier = ui->SC_fc_SpinBox->value();}
    if(senderObj == ui->SC_SymRate_SpinBox) { UISource.SC_symrate = ui->SC_SymRate_SpinBox->value(); }
    if(senderObj == ui->SC_Rolloff_doubleSpinBox) { UISource.SC_rolloff = ui->SC_Rolloff_doubleSpinBox->value(); }
    if(senderObj == ui->SC_FilterLength_SpinBox) { UISource.SC_filter_length = ui->SC_FilterLength_SpinBox->value(); }
    if(senderObj == ui->SC_FilterType_ComboBox) { UISource.SC_FilterType = ui->SC_FilterType_ComboBox->currentText(); }

    if(senderObj == ui->OFDM_fc_SpinBox) { UISource.OFDM_f_carrier = ui->OFDM_fc_SpinBox->value(); CurrentRecalcNeeds.RecalcOfdmFc = true; }
    if(senderObj == ui->OFDM_Nfft_SpinBox) { UISource.OFDM_Nfft = ui->OFDM_Nfft_SpinBox->value(); CurrentRecalcNeeds.RecalcOfdmSig = true; }
    if(senderObj == ui->OFDM_GB_DC_SpinBox) { UISource.OFDM_GB_DC = ui->OFDM_GB_DC_SpinBox->value(); CurrentRecalcNeeds.RecalcOfdmFc = true;}
    if(senderObj == ui->OFDM_GBNyq_SpinBox) { UISource.OFDM_GB_Nyq = ui->OFDM_GBNyq_SpinBox->value(); CurrentRecalcNeeds.RecalcOfdmFc = true; }
    if(senderObj == ui->OFDM_CyclePref_SpinBox) { UISource.OFDM_cycle_prefix = ui->OFDM_CyclePref_SpinBox->value(); CurrentRecalcNeeds.RecalcOfdmFc = true; }

    if(senderObj == ui->FDMA_fc_SpinBox) { UISource.FDMA_f_carrier = ui->FDMA_fc_SpinBox->value(); }
    if(senderObj == ui->FDMA_SymRate_SpinBox) { UISource.FDMA_symrate = ui->FDMA_SymRate_SpinBox->value(); }
    if(senderObj == ui->FDMA_NumCarriers_SpinBox) { UISource.FDMA_num_subcarriers = ui->FDMA_NumCarriers_SpinBox->value(); }
    if(senderObj == ui->FDMA_StepCarrier_SpinBox) { UISource.FDMA_step_carrier = ui->FDMA_StepCarrier_SpinBox->value(); }

    if(senderObj == ui->Oversapmling_SpinBox) { UISource.oversampling = ui->Oversapmling_SpinBox->value(); CurrentRecalcNeeds.RecalcOfdmFc = true; }
    if(senderObj == ui->Fs_SpinBox) { UISource.fs = ui->Fs_SpinBox->value(); CurrentRecalcNeeds.RecalcOfdmFc = true; }
    if(senderObj == ui->SNRSig_SpinBox) { UISource.SNRSig = ui->SNRSig_SpinBox->value(); CurrentRecalcNeeds.RecalcNoiseSig = true;}

    MySigProc.DataUpdate(UISource);
    MakeMainCalcAndPlot();
}

void MainWindow::FirstDataUpdate()
{
    if(ui->ModTypeComboBox->currentText() == "BPSK") UISource.M = 2;
    else if(ui->ModTypeComboBox->currentText() == "QPSK") UISource.M = 4;
    else if(ui->ModTypeComboBox->currentText() == "16QAM") UISource.M = 16;
    else if(ui->ModTypeComboBox->currentText() == "64QAM") UISource.M = 64;

    UISource.ModType = ui->ModTypeComboBox->currentText();
    UISource.SNRSymdB = ui->SNRSymSpinBox->value();
    UISource.NumSym = ui->NumSymSpinBox->value();
    UISource.SigType = ui->SignalTypeComboBox->currentText();
    UISource.SC_f_carrier = ui->SC_fc_SpinBox->value();
    UISource.SC_symrate = ui->SC_SymRate_SpinBox->value();
    UISource.SC_rolloff = ui->SC_Rolloff_doubleSpinBox->value();
    UISource.SC_filter_length = ui->SC_FilterLength_SpinBox->value();
    UISource.SC_FilterType = ui->SC_FilterType_ComboBox->currentText();
    UISource.OFDM_f_carrier = ui->OFDM_fc_SpinBox->value();
    UISource.OFDM_Nfft = ui->OFDM_Nfft_SpinBox->value();
    UISource.OFDM_GB_DC = ui->OFDM_GB_DC_SpinBox->value();
    UISource.OFDM_GB_Nyq = ui->OFDM_GBNyq_SpinBox->value();
    UISource.OFDM_cycle_prefix = ui->OFDM_CyclePref_SpinBox->value();
    UISource.FDMA_f_carrier = ui->FDMA_fc_SpinBox->value();
    UISource.FDMA_symrate = ui->FDMA_SymRate_SpinBox->value();
    UISource.FDMA_num_subcarriers = ui->FDMA_NumCarriers_SpinBox->value();
    UISource.FDMA_step_carrier = ui->FDMA_StepCarrier_SpinBox->value();
    UISource.oversampling = ui->Oversapmling_SpinBox->value();
    UISource.fs = ui->Fs_SpinBox->value();
    UISource.SNRSig = ui->SNRSig_SpinBox->value();

    CurrentRecalcNeeds.init();

    MySigProc.DataUpdate(UISource);
    MakeMainCalcAndPlot();
}

void MainWindow::onPipelineItemChanged(QTreeWidgetItem* current, QTreeWidgetItem*)
{
    if (!current) return;

    QString name = current->text(0);

    if (name == "Modulation") {

        ui->settingsStack->setCurrentWidget(ui->PageModulation);    
    }
    else if (name == "SignalType")
        ui->settingsStack->setCurrentWidget(ui->PageSignalType);
    else if (name == "Power Amplifier") {
        PaCurvePlot();
        ui->settingsStack->setCurrentWidget(ui->PagePA);   
    }
}

void MainWindow::SignalTypeComboBoxTextChanged(const QString& comboxString)
{
    if (comboxString == "SC")
        ui->PagesOfSigParametersStack->setCurrentWidget(ui->SC_page);
    else if (comboxString == "OFDM")
        ui->PagesOfSigParametersStack->setCurrentWidget(ui->OFDM_page);
    else if (comboxString == "FDMA")
        ui->PagesOfSigParametersStack->setCurrentWidget(ui->FDMA_page);
}

void MainWindow::PaModelTypeComboBoxTextChanged(const QString& comboxString)
{
    if (comboxString == "Saleh")
        ui->PaSettings_StackedWidget->setCurrentWidget(ui->Saleh_model_settings);
    else if (comboxString == "Rapp")
        ui->PaSettings_StackedWidget->setCurrentWidget(ui->Rapp_model_settings);
    else if (comboxString == "Ghorbani")
        ui->PaSettings_StackedWidget->setCurrentWidget(ui->Ghorbani_model_settings);
    else if (comboxString == "Memory Polynomial")
        ui->PaSettings_StackedWidget->setCurrentWidget(ui->MP_moodel_settings);
}

//SETUP WINDOW///////////////////////////////////////////
/// MainWindow::setupAdaptiveWindow
////////////////////////////////////////////////////////
void MainWindow::setupAdaptiveWindow()
{
    // Получаем текущий экран
    QScreen *screen = QGuiApplication::primaryScreen();

    // Берем 3/4 от доступной области экрана
    QRect screenGeometry = screen->availableGeometry();
    int width = screenGeometry.width() * 3 / 4;
    int height = screenGeometry.height() * 3 / 4;

    // Центрируем окно
    int x = screenGeometry.x() + (screenGeometry.width() - width) / 2;
    int y = screenGeometry.y() + (screenGeometry.height() - height) / 2;

    // Устанавливаем размер и позицию
    setGeometry(x, y, width, height);

    // Позволяем окну менять размер
    setMinimumSize(800, 600); // минимальный размер
}

void MainWindow::setupSplitter()
{
    // Устанавливаем пропорции разделителя: 30% дерево, 70% остальное
    QList<int> sizes;
    sizes << width() * 0.15 << width() * 0.35 << width() * 0.5;
    ui->splitter->setSizes(sizes);

    ui->splitter->setStretchFactor(0, 1);  // pipelineTree
    ui->splitter->setStretchFactor(1, 2);  // settingsStack
    ui->splitter->setStretchFactor(2, 3);  // plotArea
}

void MainWindow::setupLabels()
{
    ui->SalehFormulaLabel->setText(
        "A(r) = αₐ·r / (1 + βₐ·r²)\n\n"
        "Φ(r) = αᵩ·r² / (1 + βᵩ·r²)\n"
        );
    ui->SalehCoef1_label->setText(
        "<span style='font-size:14pt;'>"
        "&alpha;<sub>a</sub>"
        "</span>");
    ui->SalehCoef2_label->setText(
        "<span style='font-size:14pt;'>"
        "&beta;<sub>a</sub>"
        "</span>");
    ui->SalehCoef3_label->setText(
        "<span style='font-size:14pt;'>"
        "&alpha;<sub>&phi;</sub>"
        "</span>");
    ui->SalehCoef4_label->setText(
        "<span style='font-size:14pt;'>"
        "&beta;<sub>&phi;</sub>"
        "</span>");

    ui->RappAsatCoef_label->setText(
        "<span style='font-size:14pt;'>"
        "Точка насыщения A<sub>sat;</sub>"
        "</span>");

    ui->RappFormulaLabel->setText(
        "k·|A|·singA / (1 + (k|A| / s)²ᴾ)⅟²ᴾ\n\n"
        "g(A) = 0\n");

    ui->GhorbaniCoef1_label->setText(
        "<span style='font-size:14pt;'>"
        "&alpha;<sub>a;</sub>"
        "</span>");
    ui->GhorbaniCoef2_label->setText(
        "<span style='font-size:14pt;'>"
        "&beta;<sub>a;</sub>"
        "</span>");
    ui->GhorbaniCoef3_label->setText(
        "<span style='font-size:14pt;'>"
        "&gamma;<sub>a;</sub>"
        "</span>");
    ui->GhorbaniCoef4_label->setText(
        "<span style='font-size:14pt;'>"
        "&alpha;<sub>p;</sub>"
        "</span>");
    ui->GhorbaniCoef5_label->setText(
        "<span style='font-size:14pt;'>"
        "&beta;<sub>p;</sub>"
        "</span>");
    ui->GhorbaniCoef6_label->setText(
        "<span style='font-size:14pt;'>"
        "&gamma;<sub>p;</sub>"
        "</span>");
}

void MainWindow::setupMainPipelineTree()
{
    ui->pipelineTree->expandAll();
    ui->pipelineTree->topLevelItem(0)->child(0)->setSelected(true);
    //Connecting pipeline tree with settings and setting starting page
    connect(ui->pipelineTree, &QTreeWidget::currentItemChanged, this, &MainWindow::onPipelineItemChanged);
    ui->settingsStack->setCurrentWidget(ui->PageModulation);
}

void MainWindow::SetupSelectedGraphsListWidget()
{
    ui->SelectedGraphs_ListWidget->item(0)->setSelected(true);
    connect(ui->SelectedGraphs_ListWidget, &QListWidget::currentItemChanged, this, &MainWindow::onGraphsListItemChanged);
    ui->GraphsListstackedWidget->setCurrentWidget(ui->ConstellationGraphsPage);
}

void MainWindow::SetupMainLogicWork()
{
    connect(ui->ModTypeComboBox, &QComboBox::currentTextChanged, this, &MainWindow::DataUpdate);
    connect(ui->NumSymSpinBox, &QSpinBox::valueChanged, this, &MainWindow::DataUpdate);
    connect(ui->SNRSymSpinBox, &QSpinBox::valueChanged, this, &MainWindow::DataUpdate);
    connect(ui->SC_fc_SpinBox, &QSpinBox::valueChanged, this, &MainWindow::DataUpdate);
    connect(ui->SC_SymRate_SpinBox, &QSpinBox::valueChanged, this, &MainWindow::DataUpdate);
    connect(ui->SC_Rolloff_doubleSpinBox, &QDoubleSpinBox::valueChanged, this, &MainWindow::DataUpdate);
    connect(ui->SC_FilterLength_SpinBox, &QSpinBox::valueChanged, this, &MainWindow::DataUpdate);
    connect(ui->SC_FilterType_ComboBox, &QComboBox::currentTextChanged, this, &MainWindow::DataUpdate);
    connect(ui->OFDM_fc_SpinBox, &QSpinBox::valueChanged, this, &MainWindow::DataUpdate);
    connect(ui->OFDM_Nfft_SpinBox, &QSpinBox::valueChanged, this, &MainWindow::DataUpdate);
    connect(ui->OFDM_GB_DC_SpinBox, &QSpinBox::valueChanged, this, &MainWindow::DataUpdate);
    connect(ui->OFDM_GBNyq_SpinBox, &QSpinBox::valueChanged, this, &MainWindow::DataUpdate);
    connect(ui->FDMA_fc_SpinBox, &QSpinBox::valueChanged, this, &MainWindow::DataUpdate);
    connect(ui->FDMA_SymRate_SpinBox, &QSpinBox::valueChanged, this, &MainWindow::DataUpdate);
    connect(ui->FDMA_NumCarriers_SpinBox, &QSpinBox::valueChanged, this, &MainWindow::DataUpdate);
    connect(ui->FDMA_StepCarrier_SpinBox, &QSpinBox::valueChanged, this, &MainWindow::DataUpdate);
    connect(ui->SignalTypeComboBox, &QComboBox::currentIndexChanged, this, &MainWindow::DataUpdate);
    connect(ui->Oversapmling_SpinBox, &QSpinBox::valueChanged, this, &MainWindow::DataUpdate);
    connect(ui->Fs_SpinBox, &QSpinBox::valueChanged, this, &MainWindow::DataUpdate);
    connect(ui->SNRSig_SpinBox, &QSpinBox::valueChanged, this, &MainWindow::DataUpdate);
}

void MainWindow::onGraphsListItemChanged(QListWidgetItem* current, QListWidgetItem*) {
    if (!current) return;

    QString name = current->text();

    if (name == "    Constellation    ") {

        ui->GraphsListstackedWidget->setCurrentWidget(ui->ConstellationGraphsPage);
    }
    else if (name == "    PSD    ")
        ui->GraphsListstackedWidget->setCurrentWidget(ui->PSDGraphsPage);
    else if (name == "    DPD learning    ") {
        PaCurvePlot();
        ui->GraphsListstackedWidget->setCurrentWidget(ui->DPDlearningPage);
    }
    else if (name == "    Custom    ") {
        PaCurvePlot();
        ui->GraphsListstackedWidget->setCurrentWidget(ui->CustomGraphsPage);
    }
    else if (name == "    TimeDomain    ") {
        PaCurvePlot();
        ui->GraphsListstackedWidget->setCurrentWidget(ui->TimeDomainPage);
    }
}

void MainWindow::setupPACurvePlotting() {
    //Connecting chosen PA model type with settings of this model and setting starting page
    connect(ui->PAModel_ComboBox, &QComboBox::currentTextChanged, this, &MainWindow::PaModelTypeComboBoxTextChanged);
    ui->PaSettings_StackedWidget->setCurrentWidget(ui->Saleh_model_settings);

    connect(ui->PAModel_ComboBox, &QComboBox::currentTextChanged, this, &MainWindow::PaCurvePlot);
    connect(ui->IBO_SpinBox, &QSpinBox::valueChanged, this, &MainWindow::PaCurvePlot);
    connect(ui->LinearGain_SpinBox, &QSpinBox::valueChanged, this, &MainWindow::PaCurvePlot);

    connect(ui->SalehCoef1_doubleSpinBox, &QDoubleSpinBox::valueChanged, this, &MainWindow::PaCurvePlot);
    connect(ui->SalehCoef2_doubleSpinBox, &QDoubleSpinBox::valueChanged, this, &MainWindow::PaCurvePlot);
    connect(ui->SalehCoef3_doubleSpinBox, &QDoubleSpinBox::valueChanged, this, &MainWindow::PaCurvePlot);
    connect(ui->SalehCoef4_doubleSpinBox, &QDoubleSpinBox::valueChanged, this, &MainWindow::PaCurvePlot);

    connect(ui->RappAsatCoef_doubleSpinBox, &QDoubleSpinBox::valueChanged, this, &MainWindow::PaCurvePlot);
    connect(ui->RappPCoef_doubleSpinBox, &QDoubleSpinBox::valueChanged, this, &MainWindow::PaCurvePlot);

    connect(ui->GhorbaniCoef1_doubleSpinBox, &QDoubleSpinBox::valueChanged, this, &MainWindow::PaCurvePlot);
    connect(ui->GhorbaniCoef2_doubleSpinBox, &QDoubleSpinBox::valueChanged, this, &MainWindow::PaCurvePlot);
    connect(ui->GhorbaniCoef3_doubleSpinBox, &QDoubleSpinBox::valueChanged, this, &MainWindow::PaCurvePlot);
    connect(ui->GhorbaniCoef4_doubleSpinBox, &QDoubleSpinBox::valueChanged, this, &MainWindow::PaCurvePlot);
    connect(ui->GhorbaniCoef5_doubleSpinBox, &QDoubleSpinBox::valueChanged, this, &MainWindow::PaCurvePlot);
    connect(ui->GhorbaniCoef6_doubleSpinBox, &QDoubleSpinBox::valueChanged, this, &MainWindow::PaCurvePlot);
}

void MainWindow::setupPaCurveToolButton() {
    // Create menu with selectable options
    QMenu *menu = new QMenu(ui->PACurveToolButton);
    QFont font("Times New Roman", 14);

    // Add actions with checkable items
    QAction *action1 = menu->addAction("Linear");
    action1->setFont(font);
    QAction *action2 = menu->addAction("dB");
    action2->setFont(font);
    QAction *action3 = menu->addAction("Absolute");
    action3->setFont(font);
    QAction *action4 = menu->addAction("Normalized");
    action4->setFont(font);

    // Make actions checkable
    action1->setCheckable(true);
    action2->setCheckable(true);
    action3->setCheckable(true);
    action4->setCheckable(true);

    // Create an action group for exclusive selection
    QActionGroup *actionGroup1 = new QActionGroup(menu);
    QActionGroup *actionGroup2 = new QActionGroup(menu);
    actionGroup1->addAction(action1);
    actionGroup1->addAction(action2);
    actionGroup1->setExclusive(true); // Only one can be selected

    actionGroup2->addAction(action3);
    actionGroup2->addAction(action4);
    actionGroup2->setExclusive(true); // Only one can be selected

    // Set default selection
    action1->setChecked(true);
    action4->setChecked(true);

    // Set menu to button
    ui->PACurveToolButton->setMenu(menu);
    ui->PACurveToolButton->setPopupMode(QToolButton::MenuButtonPopup); // Shows arrow for menu

    connect(action1, &QAction::triggered, this, &MainWindow::PaCurvePlot);
    connect(action2, &QAction::triggered, this, &MainWindow::PaCurvePlot);
    connect(action3, &QAction::triggered, this, &MainWindow::PaCurvePlot);
    connect(action4, &QAction::triggered, this, &MainWindow::PaCurvePlot);
}
