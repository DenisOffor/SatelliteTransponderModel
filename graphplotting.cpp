    #include "graphplotting.h"

GraphPlotting::GraphPlotting(SignalProcessing* sig_proc, Ui::MainWindow *ui, QObject *parent)
    : QObject(parent)
    , Local_ui_copy(ui)
    , Local_copy_SigProc(sig_proc){
    PACurveType = "Static";
}

GraphPlotting::~GraphPlotting() {
    delete plotIdealSymConst;
    delete plotPaCurve;
    for(int i = 0; i < 6; i++) {
        delete plotsOfConstellations[i];
        delete plotsOfTimeDomain[i];
        delete plotsOfPSD[i];
        delete plotsOfDPD[i];
    }
}

void GraphPlotting::PaCurvePlot()
{
    plotPaCurve->graph(0)->data()->clear();
    plotPaCurve->graph(1)->data()->clear();
    plotPaCurve->graph(2)->data()->clear();
    plotPaCurve->graph(3)->data()->clear();

    if(PACurveType == "Static")
        PlotStaticPaCurve(Local_copy_SigProc->getPaCurve(), Local_ui_copy->PACurveToolButton->menu()->actions());
    else if(PACurveType == "Scatter")
        PlotScatterPaCurve(Local_copy_SigProc->getTimeSignal());
}

void GraphPlotting::init(std::vector<QWidget*> SetupGraphWidgets, std::vector<QWidget*> ConstellationsGraphWidgets,
                         std::vector<QWidget*> TimeDomainGraphWidgets, std::vector<QWidget*> PSDGraphWidgets, std::vector<QWidget*> DPDLearnGraphWidgets) {
    InitializeIdealSymConstPlot(SetupGraphWidgets[0]);
    InitializePaCurvePlot(SetupGraphWidgets[1]);
    InitializeConstellationsPlotting(ConstellationsGraphWidgets);
    InitializeTimeDomainPlotting(TimeDomainGraphWidgets);
    InitializePSDPlotting(PSDGraphWidgets);
    InitializeDPDLearnPlotting(DPDLearnGraphWidgets);
    setupGraphActions();
}

void GraphPlotting::PlotIdealSymConstellation(const QString ModType) {
    QVector<double> x, y;

    if (ModType == "BPSK") {
        // Две точки: (-1, 0) и (1, 0)
        x << -1.0 << 1.0;
        y << 0.0 << 0.0;
    }
    else if (ModType == "QPSK") {
        // Четыре точки
        double val = 0.707; // 1/√2
        x << -val << -val << val << val;
        y << -val << val << -val << val;
    }
    else if (ModType == "16QAM") {
        // 16 точек в сетке 4x4
        QVector<double> levels = {-1.0, -0.333, 0.333, 1.0};
        for (double i : levels) {
            for (double q : levels) {
                x << i;
                y << q;
            }
        }
    }
    else if (ModType == "64QAM") {
        // Сетка 8x8
        QVector<double> levels;
        for (int i = -7; i <= 7; i += 2) {
            levels << i / 7.0; // Нормализуем к [-1, 1]
        }
        for (double i : levels) {
            for (double q : std::as_const(levels)) {
                x << i;
                y << q;
            }
        }
    }

    // Устанавливаем данные
    plotIdealSymConst->graph(0)->setData(x, y);

    // Перерисовываем
    plotIdealSymConst->replot();
}

void GraphPlotting::PlotConstellationsPlots(const Symbols& MySymbols)
{
    QVector<double> x, y;
    if(!MySymbols.tr_sym_clean.empty()) {
        for (const auto& symbol : MySymbols.tr_sym_clean) {
            x.append(symbol.real());  // X = действительная часть
            y.append(symbol.imag());  // Y = мнимая часть
        }
        plotsOfConstellations[0]->graph(0)->setData(x, y);
        plotsOfConstellations[0]->replot();
    }

    if(!MySymbols.tr_sym_noisy.empty()) {
        x.clear(); y.clear();
        for (const auto& symbol : MySymbols.tr_sym_noisy) {
            x.append(symbol.real());  // X = действительная часть
            y.append(symbol.imag());  // Y = мнимая часть
        }
        plotsOfConstellations[1]->graph(0)->setData(x, y);
        plotsOfConstellations[1]->replot();
    }

    if(!MySymbols.rec_sym_noisy_with_DPD.empty()) {
        x.clear(); y.clear();
        for (const auto& symbol : MySymbols.rec_sym_noisy_with_DPD) {
            x.append(symbol.real());  // X = действительная часть
            y.append(symbol.imag());  // Y = мнимая часть
        }
        plotsOfConstellations[3]->graph(0)->setData(x, y);
        plotsOfConstellations[3]->replot();
    }

    if(!MySymbols.rec_sym_noisy.empty()) {
        x.clear(); y.clear();
        for (const auto& symbol : MySymbols.rec_sym_noisy) {
            x.append(symbol.real());  // X = действительная часть
            y.append(symbol.imag());  // Y = мнимая часть
        }
        plotsOfConstellations[4]->graph(0)->setData(x, y);
        plotsOfConstellations[4]->replot();
    }
}

void GraphPlotting::PlotTimeDomainPlots(const GlobalResults &CurrentRes, bool rescale)
{
    if(!CurrentRes.tx_sig.empty()) {
        QVector<double> y_tx, y_pa;
        for (const auto& symbol : CurrentRes.tx_sig) y_tx.append(symbol.real());
        for (const auto& symbol : CurrentRes.pa_sig) y_pa.append(symbol.real());

        QVector<double> scaledX;
        double scaleFactor = 1e6;  // 10^6
        scaledX.reserve(CurrentRes.time.size());
        for (double x : CurrentRes.time)
            scaledX.append(x * scaleFactor);

        QCPRange ranges[6];
        for(int i = 0; i < 6; ++i)
            ranges[i] = plotsOfTimeDomain[i]->xAxis->range();

        plotsOfTimeDomain[0]->graph(0)->setData(scaledX, y_tx);
        plotsOfTimeDomain[1]->graph(0)->setData(scaledX, y_pa);

        if(rescale) {
            for(int i = 0; i < 6; ++i)
                plotsOfTimeDomain[i]->graph(0)->rescaleAxes();
        }
        else {
            for(int i = 0; i < 6; ++i)
                plotsOfTimeDomain[i]->xAxis->setRange(ranges[i]);
        }

        for(int i = 0; i < 6; ++i) plotsOfTimeDomain[i]->replot();
    }
}

void GraphPlotting::PlotPSDPlots(const std::vector<std::vector<double>>& PSDs, const std::vector<std::vector<double>>& freqs)
{
    QVector<QVector<double>> qPSDs;
    QVector<QVector<double>> qfreqs;
    for (const auto& psd : PSDs) { qPSDs.append(QVector<double>(psd.begin(), psd.end())); }
    for (const auto& freqs : freqs) { qfreqs.append(QVector<double>(freqs.begin(), freqs.end())); }

    if(!PSDs.empty() && !freqs.empty()) {
        plotsOfPSD[0]->graph(0)->setData(qfreqs[0], qPSDs[0]);
        plotsOfPSD[1]->graph(0)->setData(qfreqs[0], qPSDs[1]);
        plotsOfPSD[0]->graph(0)->rescaleAxes();
        plotsOfPSD[1]->graph(0)->rescaleAxes();
        plotsOfPSD[0]->replot();
        plotsOfPSD[1]->replot();

        plotsOfPSD[3]->graph(0)->setData(qfreqs[0], qPSDs[0]);
        plotsOfPSD[3]->graph(1)->setData(qfreqs[0], qPSDs[1]);
        plotsOfPSD[3]->graph(0)->rescaleAxes();
        plotsOfPSD[3]->replot();


        plotsOfPSD[4]->graph(0)->setData(qfreqs[0], qPSDs[2]);
        plotsOfPSD[4]->graph(1)->setData(qfreqs[0], qPSDs[3]);
        plotsOfPSD[4]->graph(0)->rescaleAxes();
        plotsOfPSD[4]->replot();
    }
}

QString GraphPlotting::GetPaCurveType()
{
    return PACurveType;
}

void GraphPlotting::PlotStaticPaCurve(PaCurve& PACurve, const QList<QAction*> actions)
{
    QVector<double> x, y;
    QVector<double> wp_x, wp_y;
    QString xLabel, yLabel;
    x.resize(PACurve.point_num);
    y.resize(PACurve.point_num);

    if (actions[0]->isChecked() && actions[3]->isChecked()) {
        x = QVector<double>((PACurve.P_in_norm.linear.begin() + PA_CURVE_BIAS),
                            PACurve.P_in_norm.linear.end());
        y = QVector<double>((PACurve.P_out_norm.linear.begin() + PA_CURVE_BIAS),
                            PACurve.P_out_norm.linear.end());
        wp_y  = QVector<double>(PACurve.Working_point_linear_norm.y.begin(),
                            PACurve.Working_point_linear_norm.y.end());
        wp_x  = QVector<double>(PACurve.Working_point_linear_norm.x.begin(),
                            PACurve.Working_point_linear_norm.x.end());
        xLabel = "Pвх/Pнас, Вт";
        yLabel = "Pвх/Pнас, Вт";
    }
    else if (actions[0]->isChecked() && actions[2]->isChecked()) {
        x = QVector<double>((PACurve.P_in_abs.linear.begin() + PA_CURVE_BIAS),
                            PACurve.P_in_abs.linear.end());
        y = QVector<double>((PACurve.P_out_abs.linear.begin() + PA_CURVE_BIAS),
                            PACurve.P_out_abs.linear.end());
        wp_y  = QVector<double>(PACurve.Working_point_linear_abs.y.begin(),
                               PACurve.Working_point_linear_abs.y.end());
        wp_x  = QVector<double>(PACurve.Working_point_linear_abs.x.begin(),
                               PACurve.Working_point_linear_abs.x.end());
        xLabel = "Pвх, Вт";
        yLabel = "Pвых, Вт";
    }
    else if (actions[1]->isChecked() && actions[3]->isChecked()) {
        x = QVector<double>((PACurve.P_in_norm.dB.begin() + PA_CURVE_BIAS),
                            PACurve.P_in_norm.dB.end());
        y = QVector<double>((PACurve.P_out_norm.dB.begin() + PA_CURVE_BIAS),
                            PACurve.P_out_norm.dB.end());
        wp_y  = QVector<double>(PACurve.Working_point_dB_norm.y.begin(),
                               PACurve.Working_point_dB_norm.y.end());
        wp_x  = QVector<double>(PACurve.Working_point_dB_norm.x.begin(),
                               PACurve.Working_point_dB_norm.x.end());
        xLabel = "Pвх/Pнас, дБВт";
        yLabel = "Pвых/Pнас, дБВт";
    }
    else {
        x = QVector<double>((PACurve.P_in_abs.dB.begin() + PA_CURVE_BIAS),
                            PACurve.P_in_abs.dB.end());
        y = QVector<double>((PACurve.P_out_abs.dB.begin() + PA_CURVE_BIAS),
                            PACurve.P_out_abs.dB.end());
        wp_y  = QVector<double>(PACurve.Working_point_dB_abs.y.begin(),
                               PACurve.Working_point_dB_abs.y.end());
        wp_x  = QVector<double>(PACurve.Working_point_dB_abs.x.begin(),
                               PACurve.Working_point_dB_abs.x.end());
        xLabel = "Pвх, дБВт";
        yLabel = "Pвых, дБВт";
    }

    QVector<double> Phi = QVector<double>((PACurve.Phi.begin() + PA_CURVE_BIAS), PACurve.Phi.end());
    QVector<double> Phi_work = QVector<double>(PACurve.Phi_work_grad.begin(), PACurve.Phi_work_grad.end());
    plotPaCurve->graph(0)->setData(x, y);
    plotPaCurve->graph(1)->setData(x, Phi);
    plotPaCurve->graph(2)->setData(wp_x, wp_y);
    plotPaCurve->graph(3)->setData(wp_x, Phi_work);

    plotPaCurve->xAxis->setLabel(xLabel);
    plotPaCurve->yAxis->setLabel(yLabel);

    plotPaCurve->graph(0)->rescaleAxes();
    plotPaCurve->graph(1)->rescaleAxes();

    auto yRange = plotPaCurve->yAxis->range();
    double padding = (yRange.upper - yRange.lower) * 0.05;
    plotPaCurve->yAxis->setRange(yRange.lower, yRange.upper + padding);

    plotPaCurve->replot();
}

void GraphPlotting::PlotScatterPaCurve(GlobalResults &res)
{
    QVector<double> x, y;
    QVector<double> Phi;

    for(int i = 0; i < res.tx_sig.size(); i += 10) {
        x.append(std::abs(res.tx_sig[i]));
        y.append(std::abs(res.pa_sig[i]));
        Phi.append(std::arg(res.pa_sig[i]) - std::arg(res.tx_sig[i]) * 180 / 3.14);
    }

    plotPaCurve->graph(2)->setData(x, y);
    plotPaCurve->graph(3)->setData(x, Phi);

    plotPaCurve->graph(2)->rescaleAxes();
    plotPaCurve->yAxis2->setRange(0, 45);
    auto yRange = plotPaCurve->yAxis->range();
    double padding = (yRange.upper - yRange.lower) * 0.05;
    plotPaCurve->yAxis->setRange(yRange.lower, yRange.upper + padding);

    plotPaCurve->replot();
}

void GraphPlotting::PlotScatterDPDLearn(GlobalResults &res)
{
    QVector<double> x, y;
    QVector<double> Phi;

    for(int i = 0; i < res.tx_sig.size(); i += 3) {
        x.append(std::abs(res.tx_sig[i]));
        y.append(std::abs(res.pa_sig[i]));
        //Phi.append(std::arg(res.pa_sig[i]) - std::arg(res.tx_sig[i]) * 180 / 3.14);
    }
    plotsOfDPD[1]->graph(0)->setData(x, y);
    //plotsOfDPD[1]->graph(1)->setData(x, Phi);
    plotsOfDPD[1]->graph(0)->rescaleAxes();
    //plotsOfDPD[1]->yAxis2->setRange(0, 45);
    auto yRange = plotsOfDPD[1]->yAxis->range();
    double padding = (yRange.upper - yRange.lower) * 0.05;
    plotsOfDPD[1]->yAxis->setRange(yRange.lower, yRange.upper + padding);
    plotsOfDPD[1]->replot();

    x.clear(); y.clear(); Phi.clear();
    for(int i = 0; i < res.tx_sig.size(); i += 3) {
        x.append(std::abs(res.tx_sig[i]));
        y.append(std::abs(res.tx_sig[i]));
        //Phi.append(std::arg(res.tx_plus_dpd_sig[i]) - std::arg(res.tx_sig[i]) * 180 / 3.14);
    }
    plotsOfDPD[0]->graph(0)->setData(x, y);
    //plotsOfDPD[3]->graph(1)->setData(x, Phi);
    plotsOfDPD[0]->graph(0)->rescaleAxes();
    //plotsOfDPD[3]->yAxis2->setRange(0, 45);
    yRange = plotsOfDPD[0]->yAxis->range();
    padding = (yRange.upper - yRange.lower) * 0.05;
    plotsOfDPD[0]->yAxis->setRange(yRange.lower, yRange.upper + padding);
    plotsOfDPD[0]->replot();

    x.clear(); y.clear(); Phi.clear();
    for(int i = 0; i < res.tx_sig.size(); i += 3) {
        x.append(std::abs(res.tx_sig[i]));
        y.append(std::abs(res.tx_plus_dpd_sig[i]));
        //Phi.append(std::arg(res.tx_plus_dpd_sig[i]) - std::arg(res.tx_sig[i]) * 180 / 3.14);
    }
    plotsOfDPD[3]->graph(0)->setData(x, y);
    //plotsOfDPD[3]->graph(1)->setData(x, Phi);
    plotsOfDPD[3]->graph(0)->rescaleAxes();
    //plotsOfDPD[3]->yAxis2->setRange(0, 45);
    yRange = plotsOfDPD[3]->yAxis->range();
    padding = (yRange.upper - yRange.lower) * 0.05;
    plotsOfDPD[3]->yAxis->setRange(yRange.lower, yRange.upper + padding);
    plotsOfDPD[3]->replot();


    x.clear(); y.clear(); Phi.clear();
    for(int i = 0; i < res.tx_sig.size(); i += 3) {
        x.append(std::abs(res.tx_sig[i]));
        y.append(std::abs(res.pa_plus_dpd_sig[i]));
        //Phi.append(std::arg(res.pa_plus_dpd_sig[i]) - std::arg(res.tx_sig[i]) * 180 / 3.14);
    }
    plotsOfDPD[4]->graph(0)->setData(x, y);
    //plotsOfDPD[4]->graph(1)->setData(x, Phi);
    plotsOfDPD[4]->graph(0)->rescaleAxes();
    //plotsOfDPD[4]->yAxis2->setRange(0, 45);
    yRange = plotsOfDPD[4]->yAxis->range();
    padding = (yRange.upper - yRange.lower) * 0.05;
    plotsOfDPD[4]->yAxis->setRange(yRange.lower, yRange.upper + padding);
    plotsOfDPD[4]->replot();
}

void GraphPlotting::InitializePaCurvePlot(QWidget* GraphWidget) {
    plotPaCurve = new QCustomPlot(GraphWidget);
    QVBoxLayout* layout = new QVBoxLayout(GraphWidget);
    layout->setContentsMargins(3, 3, 3, 3);
    plotPaCurve->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    GraphWidget->layout()->addWidget(plotPaCurve);
    plotPaCurve->addGraph(); // 0
    plotPaCurve->graph(0)->setPen(QPen(Qt::black, 3));
    plotPaCurve->addGraph(); // 1
    plotPaCurve->graph(1)->setPen(QPen(Qt::blue, 2, Qt::DashLine));
    plotPaCurve->graph(1)->setValueAxis(plotPaCurve->yAxis2);
    plotPaCurve->addGraph(); // 2
    plotPaCurve->graph(2)->setLineStyle(QCPGraph::lsNone);
    plotPaCurve->graph(2)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, Qt::black, 10));
    plotPaCurve->addGraph(); // 3
    plotPaCurve->graph(3)->setLineStyle(QCPGraph::lsNone);
    plotPaCurve->graph(3)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, Qt::blue, 10));
    plotPaCurve->graph(3)->setValueAxis(plotPaCurve->yAxis2);
    plotPaCurve->xAxis->setLabel("P_вх, Вт");
    plotPaCurve->yAxis->setLabel("P_вых, Вт");
    plotPaCurve->yAxis2->setVisible(true);
    plotPaCurve->yAxis2->setLabel("Фвых, град");
    plotPaCurve->replot();
    connect(plotPaCurve, &QCustomPlot::mousePress, this, &GraphPlotting::onPlotClick);
}

void GraphPlotting::InitializeIdealSymConstPlot(QWidget *GraphWidget)
{
    plotIdealSymConst = new QCustomPlot(GraphWidget);
    QVBoxLayout* layout = new QVBoxLayout(GraphWidget);
    layout->setContentsMargins(3, 3, 3, 3);
    plotIdealSymConst->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    GraphWidget->layout()->addWidget(plotIdealSymConst);
    plotIdealSymConst->xAxis->setRange(-1.2, 1.2);
    plotIdealSymConst->yAxis->setRange(-1.2, 1.2);
    plotIdealSymConst->addGraph();
    plotIdealSymConst->graph(0)->setLineStyle(QCPGraph::lsNone);
    plotIdealSymConst->graph(0)->setScatterStyle(QCPScatterStyle(
        QCPScatterStyle::ssCircle,
        Qt::transparent,    // Цвет заливки
        Qt::black,          // Цвет границы
        6                   // Размер
        ));
    connect(plotIdealSymConst, &QCustomPlot::mousePress, this, &GraphPlotting::onPlotClick);
    PlotIdealSymConstellation("BPSK");
}

void GraphPlotting::InitializeTimeDomainPlotting(std::vector<QWidget*> TimeDomainGraphWidgets)
{
    for(int i = 0; i < 6; i++) {
        plotsOfTimeDomain[i] = new QCustomPlot(TimeDomainGraphWidgets[i]);
        QVBoxLayout* layout = new QVBoxLayout(TimeDomainGraphWidgets[i]);
        layout->setContentsMargins(3, 3, 3, 3);
        plotsOfTimeDomain[i]->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        TimeDomainGraphWidgets[i]->layout()->addWidget(plotsOfTimeDomain[i]);
        plotsOfTimeDomain[i]->addGraph();
        plotsOfTimeDomain[i]->graph(0)->setAdaptiveSampling(true);
        plotsOfTimeDomain[i]->graph(0)->setPen(QPen(Qt::black));
        plotsOfTimeDomain[i]->graph(0)->setLineStyle(QCPGraph::lsLine);
        plotsOfTimeDomain[i]->graph(0)->setScatterStyle(QCPScatterStyle::ssNone);
        plotsOfTimeDomain[i]->xAxis->setLabel("t, мкс");
        plotsOfTimeDomain[i]->xAxis->setNumberFormat("f");
        plotsOfTimeDomain[i]->xAxis->setNumberPrecision(1);
        plotsOfTimeDomain[i]->yAxis->setLabel("Amplitude");
        plotsOfTimeDomain[i]->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
        plotsOfTimeDomain[i]->axisRect()->setRangeZoom(Qt::Horizontal);
        plotsOfTimeDomain[i]->axisRect()->setRangeDrag(Qt::Horizontal);
        connect(plotsOfTimeDomain[i], &QCustomPlot::mousePress, this, &GraphPlotting::onPlotClick);
        plotsOfTimeDomain[i]->replot();
    }
}

void GraphPlotting::InitializePSDPlotting(std::vector<QWidget*> PSDGraphWidgets)
{
    for(int i = 0; i < 6; i++) {
        plotsOfPSD[i] = new QCustomPlot(PSDGraphWidgets[i]);
        QVBoxLayout* layout = new QVBoxLayout(PSDGraphWidgets[i]);
        layout->setContentsMargins(3, 3, 3, 3);
        plotsOfPSD[i]->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        PSDGraphWidgets[i]->layout()->addWidget(plotsOfPSD[i]);
        plotsOfPSD[i]->addGraph();
        plotsOfPSD[i]->graph(0)->setAdaptiveSampling(true);
        plotsOfPSD[i]->graph(0)->setPen(QPen(Qt::black));
        plotsOfPSD[i]->graph(0)->setLineStyle(QCPGraph::lsLine);
        plotsOfPSD[i]->graph(0)->setScatterStyle(QCPScatterStyle::ssNone);
        if(i >= 3) {
            plotsOfPSD[i]->addGraph();
            if(i == 3) {
                plotsOfPSD[i]->graph(0)->setName("Перед PA");
                plotsOfPSD[i]->graph(1)->setName("После PA");
            }
            if(i == 4) {
                plotsOfPSD[i]->graph(0)->setName("PA no DPD");
                plotsOfPSD[i]->graph(1)->setName("PA with DPD");
            }
            plotsOfPSD[i]->legend->setVisible(true);
            plotsOfPSD[i]->graph(1)->setPen(QPen(Qt::red));
            plotsOfPSD[i]->graph(1)->setLineStyle(QCPGraph::lsLine);
            plotsOfPSD[i]->graph(1)->setScatterStyle(QCPScatterStyle::ssNone);
        }
        plotsOfPSD[i]->xAxis->setLabel("f, МГц");
        plotsOfPSD[i]->xAxis->setNumberFormat("f");
        plotsOfPSD[i]->xAxis->setNumberPrecision(1);
        plotsOfPSD[i]->yAxis->setLabel("Amplitude");
        plotsOfPSD[i]->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
        plotsOfPSD[i]->axisRect()->setRangeZoom(Qt::Horizontal);
        plotsOfPSD[i]->axisRect()->setRangeDrag(Qt::Horizontal);
        connect(plotsOfPSD[i], &QCustomPlot::mousePress, this, &GraphPlotting::onPlotClick);
        plotsOfPSD[i]->replot();
    }
}

void GraphPlotting::InitializeDPDLearnPlotting(std::vector<QWidget *> DPDLearnGraphWidgets)
{
    for(int i = 0; i < 6; i++) {
        plotsOfDPD[i] = new QCustomPlot(DPDLearnGraphWidgets[i]);
        QVBoxLayout* layout = new QVBoxLayout(DPDLearnGraphWidgets[i]);
        layout->setContentsMargins(3, 3, 3, 3);
        plotsOfDPD[i]->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        DPDLearnGraphWidgets[i]->layout()->addWidget(plotsOfDPD[i]);
        plotsOfDPD[i]->addGraph();
        plotsOfDPD[i]->graph(0)->setLineStyle(QCPGraph::lsNone);
        plotsOfDPD[i]->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, Qt::black, 3));
        plotsOfDPD[i]->addGraph();
        plotsOfDPD[i]->graph(1)->setLineStyle(QCPGraph::lsNone);
        plotsOfDPD[i]->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, Qt::blue, 3));
        plotsOfDPD[i]->graph(1)->setValueAxis(plotsOfDPD[i]->yAxis2);
        plotsOfDPD[i]->xAxis->setLabel("P_вх, Вт");
        plotsOfDPD[i]->yAxis->setLabel("P_вых, Вт");
        plotsOfDPD[i]->yAxis2->setVisible(true);
        plotsOfDPD[i]->yAxis2->setLabel("Фвых, град");
        plotsOfDPD[i]->replot();
        connect(plotsOfDPD[i], &QCustomPlot::mousePress, this, &GraphPlotting::onPlotClick);
    }
}

void GraphPlotting::InitializeConstellationsPlotting(std::vector<QWidget *> ConstellationsGraphWidgets)
{
    for(int i = 0; i < 6; i++) {
        plotsOfConstellations[i] = new QCustomPlot(ConstellationsGraphWidgets[i]);
        QVBoxLayout* layout = new QVBoxLayout(ConstellationsGraphWidgets[i]);
        layout->setContentsMargins(3, 3, 3, 3);
        plotsOfConstellations[i]->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        ConstellationsGraphWidgets[i]->layout()->addWidget(plotsOfConstellations[i]);
        plotsOfConstellations[i]->xAxis->setRange(-1.2, 1.2);
        plotsOfConstellations[i]->yAxis->setRange(-1.2, 1.2);
        plotsOfConstellations[i]->addGraph();
        plotsOfConstellations[i]->graph(0)->setLineStyle(QCPGraph::lsNone);
        plotsOfConstellations[i]->graph(0)->setScatterStyle(QCPScatterStyle(
            QCPScatterStyle::ssCircle,
            Qt::transparent,    // Цвет заливки
            Qt::black,          // Цвет границы
            6                   // Размер
            ));
        connect(plotsOfConstellations[i], &QCustomPlot::mousePress, this, &GraphPlotting::onPlotClick);
        plotsOfConstellations[i]->replot();
    }
}

void GraphPlotting::setupGraphActions()
{
    PaCurveGraphActions.common = IdealSymConstGraphActions.common = SymConstGraphActions.common
    = TimeDomainGraphActions.common = PSDGraphActions.common = DPDGraphActions.common =
    {
        new QAction("Сохранить как...", this),
        new QAction("Копировать", this)
    };
    //connect(PaCurveGraphActions.common[0], &QAction::triggered, this, GraphPlotting::saveGraphToFile);
    //connect(PaCurveGraphActions.common[1], &QAction::triggered, this, GraphPlotting::copyGraphToClipboard);

    PaCurveGraphActions.specific = {
        new QAction("Статические АХ и ФХ", this),
        new QAction("Scatter текущих отсчетов", this)
    };

    m_graphActions["PaCurve"] = PaCurveGraphActions;
    m_graphActions["SymConst"] = SymConstGraphActions;
    m_graphActions["IdealSymConst"] = IdealSymConstGraphActions;
    m_graphActions["TimeDomain"] = TimeDomainGraphActions;
    m_graphActions["PSD"] = PSDGraphActions;
    m_graphActions["DPD Learning"] = DPDGraphActions;
}

void GraphPlotting::updateMenuForGraphType(QMenu &menu, const QString &graphType)
{
    menu.clear();

    const auto &actions = m_graphActions[graphType];

    if (!actions.specific.isEmpty()) {
        for (QAction *action : actions.specific) {
            menu.addAction(action);
        }
        menu.addSeparator();
    }

    for (QAction *action : actions.common) {
        menu.addAction(action);
    }
}

void GraphPlotting::changeGraphType()
{
    // Диалог выбора типа графика
    QStringList types;
    types << "Линейный" << "Точечный" << "Ступенчатый" << "Сплайн";
}

void GraphPlotting::saveGraphToFile(QCustomPlot* plot)
{
    QString fileName = QFileDialog::getSaveFileName(Local_ui_copy->centralwidget,
                                                    "Сохранить график",
                                                    QDir::homePath(),
                                                    "Изображения (*.png *.jpg *.bmp *.pdf)");

    if (!fileName.isEmpty())
    {
        if (fileName.endsWith(".pdf"))
            plot->savePdf(fileName);
        else
            plot->savePng(fileName);
    }
}

void GraphPlotting::onPlotClick(QMouseEvent *event)
{
    if (event->button() == Qt::RightButton)
    {
        QObject *senderObj = sender();
        QMenu menu(Local_ui_copy->centralwidget);

        QString currentType;
        if(senderObj == plotPaCurve) currentType = "PaCurve";
        if(senderObj == plotIdealSymConst) currentType = "IdealSymConst";
        if(senderObj == plotsOfConstellations[0] || senderObj == plotsOfConstellations[1] || senderObj == plotsOfConstellations[2]
            || senderObj == plotsOfConstellations[3] || senderObj == plotsOfConstellations[4] || senderObj == plotsOfConstellations[5]) currentType = "SymConst";
        if(senderObj == plotsOfPSD[0] || senderObj == plotsOfPSD[1] || senderObj == plotsOfPSD[2]
            || senderObj == plotsOfPSD[3] || senderObj == plotsOfPSD[4] || senderObj == plotsOfPSD[5]) currentType = "PSD";
        if(senderObj == plotsOfTimeDomain[0] || senderObj == plotsOfTimeDomain[1] || senderObj == plotsOfTimeDomain[2]
            || senderObj == plotsOfTimeDomain[3] || senderObj == plotsOfTimeDomain[4] || senderObj == plotsOfTimeDomain[5]) currentType = "TimeDomain";
        if(senderObj == plotsOfDPD[0] || senderObj == plotsOfDPD[1] || senderObj == plotsOfDPD[2]
            || senderObj == plotsOfDPD[3] || senderObj == plotsOfDPD[4] || senderObj == plotsOfDPD[5]) currentType = "DPD Learning";

        updateMenuForGraphType(menu, currentType);

        QWidget *Widget = qobject_cast<QWidget*>(senderObj);
        QAction *selected = menu.exec(Widget->mapToGlobal(event->pos()));
        QCustomPlot *plotWidget = qobject_cast<QCustomPlot*>(senderObj);

        if (selected == PaCurveGraphActions.common[0])
             saveGraphToFile(plotWidget);
        else if (selected == PaCurveGraphActions.common[1])
             copyGraphToClipboard(plotWidget);

        if (selected == PaCurveGraphActions.specific[0] || selected == PaCurveGraphActions.specific[1])
            PACurveTypeChanged(selected == PaCurveGraphActions.specific[0] ? true : false);
    }
}

void GraphPlotting::PACurveTypeChanged(bool IsStatic)
{
    if(IsStatic) {
        plotPaCurve->graph(2)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, Qt::black, 10));
        plotPaCurve->graph(3)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, Qt::blue, 10));
        PACurveType = "Static";
    }
    else {
        plotPaCurve->graph(2)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, Qt::black, 3));
        plotPaCurve->graph(3)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, Qt::blue, 3));
        PACurveType = "Scatter";
    }
    if(PACurveType == "Static")
        Local_copy_SigProc->CalcPaCurve();
    PaCurvePlot();
}

void GraphPlotting::copyGraphToClipboard(QCustomPlot* plot)
{
    // Копируем график как изображение
    QClipboard *clipboard = QApplication::clipboard();

    // Вариант 1: Копировать как PNG
    QPixmap pixmap = plot->toPixmap(800, 600); // ширина, высота
    clipboard->setPixmap(pixmap);
}
