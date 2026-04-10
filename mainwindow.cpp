#include "mainwindow.h"
#include "./ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow),
    MySigProc(),
    Graphs(&MySigProc, ui, this)
{
    ui->setupUi(this);
    std::vector<QWidget*> SetupGraphWidgets;
    SetupGraphWidgets.push_back(ui->SymConstGroupBoxWidget);
    SetupGraphWidgets.push_back(ui->PAPlotCurve_Widget);

    std::vector<QWidget*> ConstellationsGraphWidgets;
    ConstellationsGraphWidgets.push_back(ui->ConstelPageGraph1_Widget);
    ConstellationsGraphWidgets.push_back(ui->ConstelPageGraph2_Widget);
    ConstellationsGraphWidgets.push_back(ui->ConstelPageGraph3_Widget);
    ConstellationsGraphWidgets.push_back(ui->ConstelPageGraph4_Widget);
    ConstellationsGraphWidgets.push_back(ui->ConstelPageGraph5_Widget);
    ConstellationsGraphWidgets.push_back(ui->ConstelPageGraph6_Widget);

    std::vector<QWidget*> TimeDomainGraphWidgets;
    TimeDomainGraphWidgets.push_back(ui->TimeDomainGraph1_Widget);
    TimeDomainGraphWidgets.push_back(ui->TimeDomainGraph2_Widget);
    TimeDomainGraphWidgets.push_back(ui->TimeDomainGraph3_Widget);
    TimeDomainGraphWidgets.push_back(ui->TimeDomainGraph4_Widget);
    TimeDomainGraphWidgets.push_back(ui->TimeDomainGraph5_Widget);
    TimeDomainGraphWidgets.push_back(ui->TimeDomainGraph6_Widget);

    std::vector<QWidget*> PSDGraphWidgets;
    PSDGraphWidgets.push_back(ui->PSDPageGraph1_Widget);
    PSDGraphWidgets.push_back(ui->PSDPageGraph2_Widget);
    PSDGraphWidgets.push_back(ui->PSDPageGraph3_Widget);
    PSDGraphWidgets.push_back(ui->PSDPageGraph4_Widget);
    PSDGraphWidgets.push_back(ui->PSDPageGraph5_Widget);
    PSDGraphWidgets.push_back(ui->PSDPageGraph6_Widget);

    std::vector<QWidget*> DPDLearningGraphWidgets;
    DPDLearningGraphWidgets.push_back(ui->DPDLearnPageGraph1_Widget);
    DPDLearningGraphWidgets.push_back(ui->DPDLearnPageGraph2_Widget);
    DPDLearningGraphWidgets.push_back(ui->DPDLearnPageGraph3_Widget);
    DPDLearningGraphWidgets.push_back(ui->DPDLearnPageGraph4_Widget);
    DPDLearningGraphWidgets.push_back(ui->DPDLearnPageGraph5_Widget);
    DPDLearningGraphWidgets.push_back(ui->DPDLearnPageGraph6_Widget);

    Graphs.init(SetupGraphWidgets, ConstellationsGraphWidgets, TimeDomainGraphWidgets, PSDGraphWidgets, DPDLearningGraphWidgets);

    //Connecting chosen signal type with settings of this signal and setting starting page
    connect(ui->SignalTypeComboBox, &QComboBox::currentTextChanged, this, &MainWindow::SignalTypeComboBoxTextChanged);
    ui->PagesOfSigParametersStack->setCurrentWidget(ui->SC_page);

    //Connecting type of modulation and his signal constellation
    connect(ui->ModTypeComboBox, &QComboBox::currentTextChanged, &Graphs, &GraphPlotting::PlotIdealSymConstellation);

    setupMainPipelineTree();
    setupLabels();
    setupAdaptiveWindow();
    setupSplitter();
    setupPaCurveToolButton();
    setupPACurvePlotting();
    SetupSelectedGraphsListWidget();
    SetupSelectedDPDType();
    SetupMainLogicWork();
    CurrentRecalcNeeds.init();
    SetupCyclePushBtn();
    SetupWorker();
    FirstDataUpdate();
}

MainWindow::~MainWindow()
{
    thread->quit();
    thread->wait();
    delete thread;
    delete worker;
    delete timer;
    delete ui;
}

void MainWindow::MakeMainCalcAndPlot()
{
    LockParChange();
    emit startSimulation(&MySigProc, CurrentRecalcNeeds);
}

void MainWindow::DataUpdate()
{
    if(ui->ModTypeComboBox->currentText() == "BPSK") UISource.M = 2;
    if(ui->ModTypeComboBox->currentText() == "QPSK") UISource.M = 4;
    if(ui->ModTypeComboBox->currentText() == "16QAM") UISource.M = 16;
    if(ui->ModTypeComboBox->currentText() == "64QAM") UISource.M = 64;


    QObject *senderObj = sender();
    CurrentRecalcNeeds.clear();

    if(senderObj == ui->DPDAutoRecalc_checkBox) {
        UISource.DPDAutoRecalc = ui->DPDAutoRecalc_checkBox->isChecked();
        if(UISource.DPDAutoRecalc) ui->DPDRecalc_pushButton->setEnabled(false);
        else ui->DPDRecalc_pushButton->setEnabled(true);
        if(UISource.DPDAutoRecalc) CurrentRecalcNeeds.DPDRecalc = true;
    }

    if(senderObj == ui->ModTypeComboBox) { UISource.ModType = ui->ModTypeComboBox->currentText(); CurrentRecalcNeeds.FullRecalc = true; }
    if(senderObj == ui->SNRSymSpinBox) { UISource.SNRSymdB = ui->SNRSymSpinBox->value(); CurrentRecalcNeeds.RecalcNoiseSym = true; }
    if(senderObj == ui->SymSNR_checkBox) { UISource.EnableSymSNR = ui->SymSNR_checkBox->isChecked(); CurrentRecalcNeeds.FullRecalc = true; }

    if(senderObj == ui->SignalTypeComboBox) { UISource.SigType = ui->SignalTypeComboBox->currentText();
        CurrentRecalcNeeds.FullRecalc = true; CurrentRecalcNeeds.TimePlotsRescale = true; if(UISource.DPDAutoRecalc) CurrentRecalcNeeds.DPDRecalc = true; }

    if(senderObj == ui->SC_fc_SpinBox) { UISource.SC_f_carrier = ui->SC_fc_SpinBox->value(); CurrentRecalcNeeds.RecalcSig = true; }
    if(senderObj == ui->SC_SymRate_SpinBox) { UISource.SC_symrate = ui->SC_SymRate_SpinBox->value(); CurrentRecalcNeeds.RecalcSig = true; }
    if(senderObj == ui->SC_Rolloff_doubleSpinBox) { UISource.SC_rolloff = ui->SC_Rolloff_doubleSpinBox->value(); CurrentRecalcNeeds.RecalcSig = true; }
    if(senderObj == ui->SC_FilterLength_SpinBox) { UISource.SC_filter_length = ui->SC_FilterLength_SpinBox->value(); CurrentRecalcNeeds.RecalcSig = true; }
    if(senderObj == ui->SC_FilterType_ComboBox) { UISource.SC_FilterType = ui->SC_FilterType_ComboBox->currentText(); CurrentRecalcNeeds.RecalcSig = true; }

    if(senderObj == ui->OFDM_fc_SpinBox) { UISource.OFDM_f_carrier = ui->OFDM_fc_SpinBox->value(); }
    if(senderObj == ui->OFDM_Nfft_SpinBox) { UISource.OFDM_Nfft = ui->OFDM_Nfft_SpinBox->value(); CurrentRecalcNeeds.RecalcSig = true; }
    if(senderObj == ui->OFDM_GB_DC_SpinBox) { UISource.OFDM_GB_DC = ui->OFDM_GB_DC_SpinBox->value(); CurrentRecalcNeeds.RecalcSig = true; }
    if(senderObj == ui->OFDM_GBNyq_SpinBox) { UISource.OFDM_GB_Nyq = ui->OFDM_GBNyq_SpinBox->value(); CurrentRecalcNeeds.RecalcSig = true; }
    if(senderObj == ui->OFDM_CyclePref_SpinBox) { UISource.OFDM_cycle_prefix = ui->OFDM_CyclePref_SpinBox->value(); CurrentRecalcNeeds.RecalcSig = true; }

    if(senderObj == ui->FDMA_fc_SpinBox) { UISource.FDMA_f_carrier = ui->FDMA_fc_SpinBox->value(); CurrentRecalcNeeds.RecalcSig = true; }
    if(senderObj == ui->FDMA_SymRate_SpinBox) { UISource.FDMA_symrate = ui->FDMA_SymRate_SpinBox->value(); CurrentRecalcNeeds.RecalcSig = true; }
    if(senderObj == ui->FDMA_NumCarriers_SpinBox) { UISource.FDMA_num_subcarriers = ui->FDMA_NumCarriers_SpinBox->value(); CurrentRecalcNeeds.RecalcSymbols = true; }
    if(senderObj == ui->FDMA_StepCarrier_SpinBox) { UISource.FDMA_step_carrier = ui->FDMA_StepCarrier_SpinBox->value(); CurrentRecalcNeeds.RecalcSig = true; }

    if(senderObj == ui->Oversapmling_SpinBox) { UISource.oversampling = ui->Oversapmling_SpinBox->value(); CurrentRecalcNeeds.RecalcSig = true; }
    if(senderObj == ui->Fs_SpinBox) { UISource.fs = ui->Fs_SpinBox->value(); CurrentRecalcNeeds.RecalcSig = true; }
    if(senderObj == ui->SNRSig_SpinBox) { UISource.SNRSig = ui->SNRSig_SpinBox->value(); CurrentRecalcNeeds.RecalcNoiseSig = true; }

    if(senderObj == ui->SalehCoef1_doubleSpinBox || senderObj == ui->SalehCoef2_doubleSpinBox ||
        senderObj == ui->SalehCoef3_doubleSpinBox || senderObj == ui->SalehCoef4_doubleSpinBox) {
        UISource.SalehCoeffs[0] = ui->SalehCoef1_doubleSpinBox->value();
        UISource.SalehCoeffs[1] = ui->SalehCoef2_doubleSpinBox->value();
        UISource.SalehCoeffs[2] = ui->SalehCoef3_doubleSpinBox->value();
        UISource.SalehCoeffs[3] = ui->SalehCoef4_doubleSpinBox->value();
        if(ui->PAModel_ComboBox->currentText() == "Saleh") {
            CurrentRecalcNeeds.PARecalc = true;
            CurrentRecalcNeeds.PaCurveReplot = true;
        }
    }
    if(senderObj == ui->RappAsatCoef_doubleSpinBox || senderObj == ui->RappPCoef_doubleSpinBox) {
        UISource.RappCoeffs[0] = ui->RappAsatCoef_doubleSpinBox->value();
        UISource.RappCoeffs[1] = ui->RappPCoef_doubleSpinBox->value();
        if(ui->PAModel_ComboBox->currentText() == "Rapp") {
            CurrentRecalcNeeds.PARecalc = true;
            CurrentRecalcNeeds.PaCurveReplot = true;
        }
    }
    if(senderObj == ui->GhorbaniCoef1_doubleSpinBox || senderObj == ui->GhorbaniCoef2_doubleSpinBox
        || senderObj == ui->GhorbaniCoef3_doubleSpinBox || senderObj == ui->GhorbaniCoef4_doubleSpinBox
               || senderObj == ui->GhorbaniCoef5_doubleSpinBox || senderObj == ui->GhorbaniCoef6_doubleSpinBox) {
        UISource.GhorbaniCoeffs[0] = ui->GhorbaniCoef1_doubleSpinBox->value();
        UISource.GhorbaniCoeffs[1] = ui->GhorbaniCoef2_doubleSpinBox->value();
        UISource.GhorbaniCoeffs[2] = ui->GhorbaniCoef3_doubleSpinBox->value();
        UISource.GhorbaniCoeffs[3] = ui->GhorbaniCoef4_doubleSpinBox->value();
        UISource.GhorbaniCoeffs[4] = ui->GhorbaniCoef5_doubleSpinBox->value();
        UISource.GhorbaniCoeffs[5] = ui->GhorbaniCoef6_doubleSpinBox->value();
        if(ui->PAModel_ComboBox->currentText() == "Ghorbani") {
            CurrentRecalcNeeds.PARecalc = true;
            CurrentRecalcNeeds.PaCurveReplot = true;
        }
    }
    if(senderObj == ui->W_FIR_alpha_doubleSpinBox) {
        UISource.W_FIRCoeffs[0] = ui->W_FIR_alpha_doubleSpinBox->value();
        if(ui->PAModel_ComboBox->currentText() == "Wiener") {
            CurrentRecalcNeeds.PARecalc = true;
            CurrentRecalcNeeds.PaCurveReplot = true;
        }
    }
    if(senderObj == ui->H_FIR_alpha_doubleSpinBox) {
        UISource.H_FIRCoeffs[0] = ui->H_FIR_alpha_doubleSpinBox->value();
        if(ui->PAModel_ComboBox->currentText() == "Hammerstein") {
            CurrentRecalcNeeds.PARecalc = true;
            CurrentRecalcNeeds.PaCurveReplot = true;
        }
    }
    if(senderObj == ui->WH_FIR_alpha1_doubleSpinBox || senderObj == ui->WH_FIR_alpha2_doubleSpinBox) {
        UISource.WH_FIRCoeffs[0] = ui->WH_FIR_alpha1_doubleSpinBox->value();
        UISource.WH_FIRCoeffs[1] = ui->WH_FIR_alpha2_doubleSpinBox->value();
        if(ui->PAModel_ComboBox->currentText() == "Wiener-Hammerstein") {
            CurrentRecalcNeeds.PARecalc = true;
            CurrentRecalcNeeds.PaCurveReplot = true;
        }
    }

    if(senderObj == ui->W_StaticNonlin_comboBox) {UISource.W_StaticNonlinModel = ui->W_StaticNonlin_comboBox->currentText();
        CurrentRecalcNeeds.PARecalc = true; CurrentRecalcNeeds.PaCurveReplot = true; }
    if(senderObj == ui->H_StaticNonlin_comboBox) {UISource.W_StaticNonlinModel = ui->W_StaticNonlin_comboBox->currentText();
        CurrentRecalcNeeds.PARecalc = true; CurrentRecalcNeeds.PaCurveReplot = true; }
    if(senderObj == ui->WH_StaticNonlin_comboBox) {UISource.W_StaticNonlinModel = ui->W_StaticNonlin_comboBox->currentText();
        CurrentRecalcNeeds.PARecalc = true; CurrentRecalcNeeds.PaCurveReplot = true; }
    if(senderObj == ui->LinearGain_SpinBox) { UISource.linear_gain_dB = ui->LinearGain_SpinBox->value();
        CurrentRecalcNeeds.PARecalc = true; CurrentRecalcNeeds.PaCurveReplot = true; CurrentRecalcNeeds.TimePlotsRescale = true;}
    if(senderObj == ui->IBO_SpinBox) { UISource.IBO_dB = ui->IBO_SpinBox->value();
        CurrentRecalcNeeds.PARecalc = true; CurrentRecalcNeeds.PaCurveReplot = true;}
    if(senderObj == ui->PAModel_ComboBox) {UISource.PAModel = ui->PAModel_ComboBox->currentText();
        CurrentRecalcNeeds.PARecalc = true; CurrentRecalcNeeds.PaCurveReplot = true; CurrentRecalcNeeds.TimePlotsRescale = true;}

    if(senderObj == ui->MP_M_spinBox) { UISource.MP_M = ui->MP_M_spinBox->value(); if(UISource.DPDAutoRecalc) CurrentRecalcNeeds.DPDRecalc = true; }
    if(senderObj == ui->MP_P_spinBox) { UISource.MP_M = ui->MP_P_spinBox->value(); if(UISource.DPDAutoRecalc) CurrentRecalcNeeds.DPDRecalc = true; }

    if(senderObj == ui->GMP_M_spinBox) { UISource.GMP_M = ui->GMP_M_spinBox->value(); if(UISource.DPDAutoRecalc) CurrentRecalcNeeds.DPDRecalc = true; }
    if(senderObj == ui->GMP_P_spinBox) { UISource.GMP_M = ui->GMP_P_spinBox->value(); if(UISource.DPDAutoRecalc) CurrentRecalcNeeds.DPDRecalc = true; }
    if(senderObj == ui->GMP_L_lag_spinBox) { UISource.GMP_L_lag = ui->GMP_L_lag_spinBox->value(); if(UISource.DPDAutoRecalc) CurrentRecalcNeeds.DPDRecalc = true; }
    if(senderObj == ui->GMP_L_lead_spinBox) { UISource.GMP_L_lead = ui->GMP_L_lead_spinBox->value(); if(UISource.DPDAutoRecalc) CurrentRecalcNeeds.DPDRecalc = true; }

    if(senderObj == ui->NormalizationType_ComboBox) { UISource.NormalizationType = ui->NormalizationType_ComboBox->currentText();
        if(UISource.DPDAutoRecalc) CurrentRecalcNeeds.DPDRecalc = true; }
    if(senderObj == ui->PredistorterType_comboBox) { UISource.PredistorterType = ui->PredistorterType_comboBox->currentText();
        if(UISource.DPDAutoRecalc) CurrentRecalcNeeds.DPDRecalc = true; }
    if(senderObj == ui->EvenP_checkBox) { UISource.Enable_even_P = ui->EvenP_checkBox->isChecked(); CurrentRecalcNeeds.DPDRecalc = true; }

    if(CurrentRecalcNeeds.PARecalc || CurrentRecalcNeeds.RecalcSig || CurrentRecalcNeeds.FullRecalc) {
        if(UISource.DPDAutoRecalc) CurrentRecalcNeeds.DPDRecalc = true;
        MySigProc.clear_OFDM_buffs();
    }

    qDebug() << "Сигнал получен от:" << senderObj->metaObject()->className();
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
    UISource.NumSym = 1000;
    UISource.EnableSymSNR = ui->SymSNR_checkBox->isChecked();
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

    UISource.SalehCoeffs.push_back(ui->SalehCoef1_doubleSpinBox->value());
    UISource.SalehCoeffs.push_back(ui->SalehCoef2_doubleSpinBox->value());
    UISource.SalehCoeffs.push_back(ui->SalehCoef3_doubleSpinBox->value());
    UISource.SalehCoeffs.push_back(ui->SalehCoef4_doubleSpinBox->value());

    UISource.RappCoeffs.push_back(ui->RappAsatCoef_doubleSpinBox->value());
    UISource.RappCoeffs.push_back(ui->RappPCoef_doubleSpinBox->value());

    UISource.GhorbaniCoeffs.push_back(ui->GhorbaniCoef1_doubleSpinBox->value());
    UISource.GhorbaniCoeffs.push_back(ui->GhorbaniCoef2_doubleSpinBox->value());
    UISource.GhorbaniCoeffs.push_back(ui->GhorbaniCoef3_doubleSpinBox->value());
    UISource.GhorbaniCoeffs.push_back(ui->GhorbaniCoef4_doubleSpinBox->value());
    UISource.GhorbaniCoeffs.push_back(ui->GhorbaniCoef5_doubleSpinBox->value());
    UISource.GhorbaniCoeffs.push_back(ui->GhorbaniCoef6_doubleSpinBox->value());

    UISource.W_FIRCoeffs.push_back(ui->W_FIR_alpha_doubleSpinBox->value());
    UISource.H_FIRCoeffs.push_back(ui->W_FIR_alpha_doubleSpinBox->value());
    UISource.WH_FIRCoeffs.push_back(ui->WH_FIR_alpha1_doubleSpinBox->value());
    UISource.WH_FIRCoeffs.push_back(ui->WH_FIR_alpha2_doubleSpinBox->value());

    UISource.linear_gain_dB = ui->LinearGain_SpinBox->value();
    UISource.IBO_dB = ui->IBO_SpinBox->value();
    UISource.PAModel = ui->PAModel_ComboBox->currentText();
    UISource.W_StaticNonlinModel = ui->W_StaticNonlin_comboBox->currentText();
    UISource.H_StaticNonlinModel = ui->W_StaticNonlin_comboBox->currentText();
    UISource.WH_StaticNonlinModel = ui->W_StaticNonlin_comboBox->currentText();

    UISource.MP_M = ui->MP_M_spinBox->value();
    UISource.MP_P = ui->MP_P_spinBox->value();

    UISource.GMP_M = ui->GMP_M_spinBox->value();
    UISource.GMP_P = ui->GMP_P_spinBox->value();
    UISource.GMP_L_lag = ui->GMP_L_lag_spinBox->value();
    UISource.GMP_L_lead = ui->GMP_L_lead_spinBox->value();

    UISource.NormalizationType = ui->NormalizationType_ComboBox->currentText();
    UISource.PredistorterType = ui->PredistorterType_comboBox->currentText();
    UISource.DPDAutoRecalc = ui->DPDAutoRecalc_checkBox->isChecked();
    UISource.Enable_even_P = ui->EvenP_checkBox->isChecked();

    CurrentRecalcNeeds.init();
    MySigProc.DataUpdate(UISource);
    emit ui->DPDRecalc_pushButton->pressed();
    MakeMainCalcAndPlot();
}

void MainWindow::SetupMainLogicWork()
{
    connect(ui->ModTypeComboBox, &QComboBox::currentTextChanged, this, &MainWindow::DataUpdate);;
    connect(ui->SNRSymSpinBox, &QSpinBox::valueChanged, this, &MainWindow::DataUpdate);
    connect(ui->SymSNR_checkBox, &QCheckBox::stateChanged, this, &MainWindow::DataUpdate);

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

    connect(ui->SalehCoef1_doubleSpinBox, &QDoubleSpinBox::valueChanged, this, &MainWindow::DataUpdate);
    connect(ui->SalehCoef2_doubleSpinBox, &QDoubleSpinBox::valueChanged, this, &MainWindow::DataUpdate);
    connect(ui->SalehCoef3_doubleSpinBox, &QDoubleSpinBox::valueChanged, this, &MainWindow::DataUpdate);
    connect(ui->SalehCoef4_doubleSpinBox, &QDoubleSpinBox::valueChanged, this, &MainWindow::DataUpdate);
    connect(ui->RappAsatCoef_doubleSpinBox, &QDoubleSpinBox::valueChanged, this, &MainWindow::DataUpdate);
    connect(ui->RappPCoef_doubleSpinBox, &QDoubleSpinBox::valueChanged, this, &MainWindow::DataUpdate);
    connect(ui->GhorbaniCoef1_doubleSpinBox, &QDoubleSpinBox::valueChanged, this, &MainWindow::DataUpdate);
    connect(ui->GhorbaniCoef2_doubleSpinBox, &QDoubleSpinBox::valueChanged, this, &MainWindow::DataUpdate);
    connect(ui->GhorbaniCoef3_doubleSpinBox, &QDoubleSpinBox::valueChanged, this, &MainWindow::DataUpdate);
    connect(ui->GhorbaniCoef4_doubleSpinBox, &QDoubleSpinBox::valueChanged, this, &MainWindow::DataUpdate);
    connect(ui->GhorbaniCoef5_doubleSpinBox, &QDoubleSpinBox::valueChanged, this, &MainWindow::DataUpdate);
    connect(ui->GhorbaniCoef6_doubleSpinBox, &QDoubleSpinBox::valueChanged, this, &MainWindow::DataUpdate);
    connect(ui->LinearGain_SpinBox, &QSpinBox::valueChanged, this, &MainWindow::DataUpdate);
    connect(ui->IBO_SpinBox, &QSpinBox::valueChanged, this, &MainWindow::DataUpdate);
    connect(ui->PAModel_ComboBox, &QComboBox::currentTextChanged, this, &MainWindow::DataUpdate);
    connect(ui->W_FIR_alpha_doubleSpinBox, &QDoubleSpinBox::valueChanged, this, &MainWindow::DataUpdate);
    connect(ui->W_StaticNonlin_comboBox, &QComboBox::currentTextChanged, this, &MainWindow::DataUpdate);
    connect(ui->H_FIR_alpha_doubleSpinBox, &QDoubleSpinBox::valueChanged, this, &MainWindow::DataUpdate);
    connect(ui->H_StaticNonlin_comboBox, &QComboBox::currentTextChanged, this, &MainWindow::DataUpdate);
    connect(ui->WH_FIR_alpha1_doubleSpinBox, &QDoubleSpinBox::valueChanged, this, &MainWindow::DataUpdate);
    connect(ui->WH_FIR_alpha2_doubleSpinBox, &QDoubleSpinBox::valueChanged, this, &MainWindow::DataUpdate);
    connect(ui->WH_StaticNonlin_comboBox, &QComboBox::currentTextChanged, this, &MainWindow::DataUpdate);

    connect(ui->MP_M_spinBox, &QSpinBox::valueChanged, this, &MainWindow::DataUpdate);
    connect(ui->MP_P_spinBox, &QSpinBox::valueChanged, this, &MainWindow::DataUpdate);

    connect(ui->GMP_M_spinBox, &QSpinBox::valueChanged, this, &MainWindow::DataUpdate);
    connect(ui->GMP_P_spinBox, &QSpinBox::valueChanged, this, &MainWindow::DataUpdate);
    connect(ui->GMP_L_lag_spinBox, &QSpinBox::valueChanged, this, &MainWindow::DataUpdate);
    connect(ui->GMP_L_lead_spinBox, &QSpinBox::valueChanged, this, &MainWindow::DataUpdate);

    connect(ui->NormalizationType_ComboBox, &QComboBox::currentTextChanged, this, &MainWindow::DataUpdate);
    connect(ui->PredistorterType_comboBox, &QComboBox::currentTextChanged, this, &MainWindow::DataUpdate);
    connect(ui->DPDRecalc_pushButton, &QPushButton::clicked, this, &MainWindow::DPDRecalcBtnClicked);  // &this->MySigProc
    connect(ui->DPDAutoRecalc_checkBox, &QCheckBox::stateChanged, this, &MainWindow::DataUpdate);
    connect(ui->EvenP_checkBox, &QCheckBox::stateChanged, this, &MainWindow::DataUpdate);
}

void MainWindow::LockParChange()
{
    ui->progressBar->setEnabled(true);
    ui->progressBar->setVisible(true);

    ui->ModTypeComboBox->setEnabled(false);
    ui->SNRSymSpinBox->setEnabled(false);
    ui->SymSNR_checkBox->setEnabled(false);

    ui->SignalTypeComboBox->setEnabled(false);
    ui->SC_fc_SpinBox->setEnabled(false);
    ui->SC_SymRate_SpinBox->setEnabled(false);
    ui->SC_Rolloff_doubleSpinBox->setEnabled(false);
    ui->SC_FilterLength_SpinBox->setEnabled(false);
    ui->SC_FilterType_ComboBox->setEnabled(false);
    ui->OFDM_fc_SpinBox->setEnabled(false);
    ui->OFDM_Nfft_SpinBox->setEnabled(false);
    ui->OFDM_GB_DC_SpinBox->setEnabled(false);
    ui->OFDM_GBNyq_SpinBox->setEnabled(false);
    ui->OFDM_CyclePref_SpinBox->setEnabled(false);
    ui->FDMA_fc_SpinBox->setEnabled(false);
    ui->FDMA_SymRate_SpinBox->setEnabled(false);
    ui->FDMA_NumCarriers_SpinBox->setEnabled(false);
    ui->FDMA_StepCarrier_SpinBox->setEnabled(false);
    ui->Oversapmling_SpinBox->setEnabled(false);
    ui->Fs_SpinBox->setEnabled(false);
    ui->SNRSig_SpinBox->setEnabled(false);

    ui->SalehCoef1_doubleSpinBox->setEnabled(false);
    ui->SalehCoef2_doubleSpinBox->setEnabled(false);
    ui->SalehCoef3_doubleSpinBox->setEnabled(false);
    ui->SalehCoef4_doubleSpinBox->setEnabled(false);

    ui->RappAsatCoef_doubleSpinBox->setEnabled(false);
    ui->RappPCoef_doubleSpinBox->setEnabled(false);

    ui->GhorbaniCoef1_doubleSpinBox->setEnabled(false);
    ui->GhorbaniCoef2_doubleSpinBox->setEnabled(false);
    ui->GhorbaniCoef3_doubleSpinBox->setEnabled(false);
    ui->GhorbaniCoef4_doubleSpinBox->setEnabled(false);
    ui->GhorbaniCoef5_doubleSpinBox->setEnabled(false);
    ui->GhorbaniCoef6_doubleSpinBox->setEnabled(false);

    ui->W_FIR_alpha_doubleSpinBox->setEnabled(false);
    ui->H_FIR_alpha_doubleSpinBox->setEnabled(false);
    ui->WH_FIR_alpha1_doubleSpinBox->setEnabled(false);
    ui->WH_FIR_alpha2_doubleSpinBox->setEnabled(false);

    ui->LinearGain_SpinBox->setEnabled(false);
    ui->IBO_SpinBox->setEnabled(false);
    ui->PAModel_ComboBox->setEnabled(false);
    ui->W_StaticNonlin_comboBox->setEnabled(false);
    ui->H_StaticNonlin_comboBox->setEnabled(false);
    ui->WH_StaticNonlin_comboBox->setEnabled(false);

    ui->MP_M_spinBox->setEnabled(false);
    ui->MP_P_spinBox->setEnabled(false);

    ui->GMP_M_spinBox->setEnabled(false);
    ui->GMP_P_spinBox->setEnabled(false);
    ui->GMP_L_lag_spinBox->setEnabled(false);
    ui->GMP_L_lead_spinBox->setEnabled(false);

    ui->NormalizationType_ComboBox->setEnabled(false);
    ui->PredistorterType_comboBox->setEnabled(false);
    ui->DPDAutoRecalc_checkBox->setEnabled(false);
    ui->EvenP_checkBox->setEnabled(false);
}

void MainWindow::UnLockParChange()
{
    ui->progressBar->setEnabled(false);
    ui->progressBar->setVisible(false);

    ui->ModTypeComboBox->setEnabled(true);
    ui->SNRSymSpinBox->setEnabled(true);
    ui->SymSNR_checkBox->setEnabled(true);

    ui->SignalTypeComboBox->setEnabled(true);
    ui->SC_fc_SpinBox->setEnabled(true);
    ui->SC_SymRate_SpinBox->setEnabled(true);
    ui->SC_Rolloff_doubleSpinBox->setEnabled(true);
    ui->SC_FilterLength_SpinBox->setEnabled(true);
    ui->SC_FilterType_ComboBox->setEnabled(true);
    ui->OFDM_fc_SpinBox->setEnabled(true);
    ui->OFDM_Nfft_SpinBox->setEnabled(true);
    ui->OFDM_GB_DC_SpinBox->setEnabled(true);
    ui->OFDM_GBNyq_SpinBox->setEnabled(true);
    ui->OFDM_CyclePref_SpinBox->setEnabled(true);
    ui->FDMA_fc_SpinBox->setEnabled(true);
    ui->FDMA_SymRate_SpinBox->setEnabled(true);
    ui->FDMA_NumCarriers_SpinBox->setEnabled(true);
    ui->FDMA_StepCarrier_SpinBox->setEnabled(true);
    ui->Oversapmling_SpinBox->setEnabled(true);
    ui->Fs_SpinBox->setEnabled(true);
    ui->SNRSig_SpinBox->setEnabled(true);

    ui->SalehCoef1_doubleSpinBox->setEnabled(true);
    ui->SalehCoef2_doubleSpinBox->setEnabled(true);
    ui->SalehCoef3_doubleSpinBox->setEnabled(true);
    ui->SalehCoef4_doubleSpinBox->setEnabled(true);

    ui->RappAsatCoef_doubleSpinBox->setEnabled(true);
    ui->RappPCoef_doubleSpinBox->setEnabled(true);

    ui->GhorbaniCoef1_doubleSpinBox->setEnabled(true);
    ui->GhorbaniCoef2_doubleSpinBox->setEnabled(true);
    ui->GhorbaniCoef3_doubleSpinBox->setEnabled(true);
    ui->GhorbaniCoef4_doubleSpinBox->setEnabled(true);
    ui->GhorbaniCoef5_doubleSpinBox->setEnabled(true);
    ui->GhorbaniCoef6_doubleSpinBox->setEnabled(true);

    ui->W_FIR_alpha_doubleSpinBox->setEnabled(true);
    ui->H_FIR_alpha_doubleSpinBox->setEnabled(true);
    ui->WH_FIR_alpha1_doubleSpinBox->setEnabled(true);
    ui->WH_FIR_alpha2_doubleSpinBox->setEnabled(true);

    ui->LinearGain_SpinBox->setEnabled(true);
    ui->IBO_SpinBox->setEnabled(true);
    ui->PAModel_ComboBox->setEnabled(true);
    ui->W_StaticNonlin_comboBox->setEnabled(true);
    ui->H_StaticNonlin_comboBox->setEnabled(true);
    ui->WH_StaticNonlin_comboBox->setEnabled(true);

    ui->MP_M_spinBox->setEnabled(true);
    ui->MP_P_spinBox->setEnabled(true);

    ui->GMP_M_spinBox->setEnabled(true);
    ui->GMP_P_spinBox->setEnabled(true);
    ui->GMP_L_lag_spinBox->setEnabled(true);
    ui->GMP_L_lead_spinBox->setEnabled(true);

    ui->NormalizationType_ComboBox->setEnabled(true);
    ui->PredistorterType_comboBox->setEnabled(true);
    ui->DPDAutoRecalc_checkBox->setEnabled(true);
    ui->EvenP_checkBox->setEnabled(true);
}

void MainWindow::handleResult()
{
    if(!timer->isActive())
        UnLockParChange();
    QElapsedTimer timer;
    timer.start();

    auto item = ui->SelectedGraphs_ListWidget->currentItem();
    QString text;
    if (item)
        text = item->text();

    if(CurrentRecalcNeeds.PaCurveReplot)
         PaCurvePlot();

    if(text == "    Constellation    ")
        Graphs.PlotConstellationsPlots(MySigProc.getSymbols());

    if(text == "    TimeDomain    ")
        Graphs.PlotTimeDomainPlots(MySigProc.getTimeSignal(), CurrentRecalcNeeds.TimePlotsRescale);

    if(text == "    PSD    ")
        Graphs.PlotPSDPlots(MySigProc.getPSDs(), MySigProc.getFreq());

    if(text == "    DPD learning    ")
        Graphs.PlotScatterDPDLearn(MySigProc.getTimeSignal());

    qDebug() << "Graphs time:" << timer.elapsed() << "ms";
    qDebug() << "\n";

    ui->BER_noDPD_doubleSpinBox->setValue(MySigProc.getTimeSignal().BER_noDPD);
    ui->BER_noDPD_label->setText(QString::number(MySigProc.getTimeSignal().BER_noDPD, 'e', 2));

    ui->BER_withDPD_doubleSpinBox->setValue(MySigProc.getTimeSignal().BER_withDPD);
    ui->BER_withDPD_label->setText(QString::number(MySigProc.getTimeSignal().BER_withDPD, 'e', 2));

    ui->EVM_noDPD_doubleSpinBox->setValue(MySigProc.getTimeSignal().EVM_noDPD);
    ui->EVM_withDPD_doubleSpinBox->setValue(MySigProc.getTimeSignal().EVM_withDPD);

    if(UISource.SigType == "OFDM")
        ui->BB_label->setText(QString::number(MySigProc.getOfdmRes().BB));
    else if(UISource.SigType == "FDMA")
        ui->BB_label->setText(QString::number(MySigProc.getFDMARes().totalBandwidth));
    else if(UISource.SigType == "SC")
        ui->BB_label->setText(QString::number(MySigProc.getSCRes().bandwidth));
}

void MainWindow::PaCurvePlot() {
    QObject *senderObj = sender();
    if (senderObj && senderObj->inherits("QAction")) {
        Graphs.PaCurvePlot();
        return;
    }

    QElapsedTimer timer;
    timer.start();
    if(ui->pipelineTree->currentItem() && ui->pipelineTree->currentItem()->text(0) == "Power Amplifier")
        if(Graphs.GetPaCurveType() == "Static")
            MySigProc.CalcPaCurve();

    Graphs.PaCurvePlot();
    qDebug() << "PaCurve time:" << timer.elapsed() << "ms";
}

void MainWindow::cycleBtnClicked()
{
    if(ui->Cycle_pushButton->isChecked()) {
        LockParChange();
        timer->start(100);
        CurrentRecalcNeeds.CycleMode = true;
        UISource.NumSym = 200;
        MySigProc.DataUpdate(UISource);
    } else {
        UnLockParChange();
        timer->stop();
        CurrentRecalcNeeds.CycleMode = false;
        CurrentRecalcNeeds.FullRecalc = true;
        UISource.NumSym = 1000;
        MySigProc.clear_OFDM_buffs();
        MySigProc.DataUpdate(UISource);
        MakeMainCalcAndPlot();
    }
}

void MainWindow::DPDRecalcBtnClicked()
{
    CurrentRecalcNeeds.DPDRecalc = true;
    MakeMainCalcAndPlot();
}

void MainWindow::CycleModeSlot()
{
    CurrentRecalcNeeds.RecalcSymbols = true;
    CurrentRecalcNeeds.DPDRecalc = false;
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
    else if (name == "Metrics") {
        ui->settingsStack->setCurrentWidget(ui->Metric_Page);
    }
    else if (name == "Predistorter") {
        ui->settingsStack->setCurrentWidget(ui->PagePredistorter);
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
    else if (comboxString == "Wiener")
        ui->PaSettings_StackedWidget->setCurrentWidget(ui->Wiener_model_settings);
    else if (comboxString == "Hammerstein")
        ui->PaSettings_StackedWidget->setCurrentWidget(ui->Hammerstein_model_settings);
    else if (comboxString == "Wiener-Hammerstein")
        ui->PaSettings_StackedWidget->setCurrentWidget(ui->WH_model_settings);
}

void MainWindow::DPDTypeComboBoxTextChange(const QString &comboxString)
{
    if (comboxString == "MP")
        ui->PredPagesType_stackedWidget->setCurrentWidget(ui->MP_Page);
    else if (comboxString == "GMP")
        ui->PredPagesType_stackedWidget->setCurrentWidget(ui->GMP_Page);
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
    sizes << width() * 0.09 << width() * 0.24 << width() * 0.67;
    ui->splitter->setSizes(sizes);

    ui->splitter->setStretchFactor(0, 1);  // pipelineTree
    ui->splitter->setStretchFactor(1, 2);  // settingsStack
    ui->splitter->setStretchFactor(2, 3);  // plotArea
}

void MainWindow::setupLabels()
{
    ui->progressBar->setEnabled(false);
    ui->progressBar->setVisible(false);
    ui->progressBar->setRange(0,0);
    ui->progressBar->setStyleSheet(
        "QProgressBar {"
        "    border: 1px solid #888;"
        "    border-radius: 5px;"
        "    background-color: #f0f0f0;"  // Светлый фон
        "}"
        "QProgressBar::chunk {"
        "    background-color: #4CAF50;"   // Зеленый цвет бегущей полоски
        "    width: 20px;"                  // Ширина сегмента (для анимации)"
        "}"
        );

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
    ui->W_FIR_Coefs_label->setText(
        "h[m]=C&alpha;<sup>m</sup>");
    ui->W_FIR_alpha_label->setText("<html>&alpha;:</html>");
}

void MainWindow::setupMainPipelineTree()
{
    ui->pipelineTree->expandAll();
    ui->pipelineTree->topLevelItem(0)->child(0)->setSelected(true);
    //Connecting pipeline tree with settings and setting starting page
    connect(ui->pipelineTree, &QTreeWidget::currentItemChanged, this, &MainWindow::onPipelineItemChanged);
    ui->settingsStack->setCurrentWidget(ui->PageModulation);

    ui->Maintain_GroupBox->setStyleSheet(
           "QGroupBox {"
           "    border: none;"           // Убираем все рамки
           "    padding-top: 0px;"       // Убираем отступ сверху для заголовка
           "    margin-top: 0px;"         // Убираем внешний отступ
           "}");
    ui->Menu_groupBox->setStyleSheet(
        "QGroupBox {"
        "    border: none;"           // Убираем все рамки
        "    padding-top: 0px;"       // Убираем отступ сверху для заголовка
        "    margin-top: 0px;"         // Убираем внешний отступ
        "}");
}

void MainWindow::SetupSelectedGraphsListWidget()
{
    ui->SelectedGraphs_ListWidget->item(0)->setSelected(true);
    connect(ui->SelectedGraphs_ListWidget, &QListWidget::currentItemChanged, this, &MainWindow::onGraphsListItemChanged);
    ui->SelectedGraphs_ListWidget->setCurrentRow(0);
    ui->GraphsListstackedWidget->setCurrentWidget(ui->ConstellationGraphsPage);
    connect(ui->SelectedGraphs_ListWidget, &QListWidget::currentItemChanged, this, &MainWindow::handleResult);
}

void MainWindow::SetupSelectedDPDType()
{
    connect(ui->PredistorterType_comboBox, &QComboBox::currentTextChanged, this, &MainWindow::DPDTypeComboBoxTextChange);
    ui->PredPagesType_stackedWidget->setCurrentWidget(ui->MP_Page);
}

void MainWindow::SetupWorker()
{
    thread = new QThread;
    worker = new Worker;

    worker->moveToThread(thread);

    connect(this, &MainWindow::startSimulation,
            worker, &Worker::process);

    connect(worker, &Worker::resultReady,
            this, &MainWindow::handleResult);

    thread->start();
}

void MainWindow::SetupCyclePushBtn()
{
    timer = new QTimer(this);
    connect(timer, &QTimer::timeout, this, &MainWindow::CycleModeSlot);

    connect(ui->Cycle_pushButton, &QPushButton::clicked, this, &MainWindow::cycleBtnClicked);
    ui->Cycle_pushButton->setCheckable(true);
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
        ui->GraphsListstackedWidget->setCurrentWidget(ui->DPDlearningPage);
    }
    else if (name == "    Custom    ") {
        ui->GraphsListstackedWidget->setCurrentWidget(ui->CustomGraphsPage);
    }
    else if (name == "    TimeDomain    ") {
        CurrentRecalcNeeds.TimePlotsRescale = true;
        ui->GraphsListstackedWidget->setCurrentWidget(ui->TimeDomainPage);
    }
}

void MainWindow::setupPACurvePlotting() {
    // //Connecting chosen PA model type with settings of this model and setting starting page
     connect(ui->PAModel_ComboBox, &QComboBox::currentTextChanged, this, &MainWindow::PaModelTypeComboBoxTextChanged);
     ui->PaSettings_StackedWidget->setCurrentWidget(ui->Saleh_model_settings);
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
