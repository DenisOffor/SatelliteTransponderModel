    #include "graphplotting.h"

GraphPlotting::GraphPlotting() {
}

GraphPlotting::~GraphPlotting() {
    delete plotIdealSymConst;
    delete plotPaCurve;
    for(int i = 0; i < 6; i++) {
        delete plotsOfConstellations[i];
        delete plotsOfTimeDomain[i];
        delete plotsOfPSD[i];
    }
}

void GraphPlotting::init(QVector<QWidget*> SetupGraphWidgets, QVector<QWidget*> ConstellationsGraphWidgets,
                         QVector<QWidget*> TimeDomainGraphWidgets, QVector<QWidget*> PSDGraphWidgets) {
    InitializeIdealSymConstPlot(SetupGraphWidgets[0]);
    InitializePaCurvePlot(SetupGraphWidgets[1]);
    InitializeConstellationsPlotting(ConstellationsGraphWidgets);
    InitializeTimeDomainPlotting(TimeDomainGraphWidgets);
    InitializePSDPlotting(PSDGraphWidgets);
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
    for (const auto& symbol : MySymbols.tr_sym_clean) {
        x.append(symbol.real());  // X = действительная часть
        y.append(symbol.imag());  // Y = мнимая часть
    }
    plotsOfConstellations[0]->graph(0)->setData(x, y);
    plotsOfConstellations[0]->replot();

    x.clear(); y.clear();
    for (const auto& symbol : MySymbols.tr_sym_noisy) {
        x.append(symbol.real());  // X = действительная часть
        y.append(symbol.imag());  // Y = мнимая часть
    }
    plotsOfConstellations[1]->graph(0)->setData(x, y);
    plotsOfConstellations[1]->replot();

    x.clear(); y.clear();
    for (const auto& symbol : MySymbols.rec_sym_noisy) {
        x.append(symbol.real());  // X = действительная часть
        y.append(symbol.imag());  // Y = мнимая часть
    }
    plotsOfConstellations[4]->graph(0)->setData(x, y);
    plotsOfConstellations[4]->replot();
}

void GraphPlotting::PlotTimeDomainPlots(const OfdmResult &CurrentOfdm)
{
    static bool first_in = true;
    if(!CurrentOfdm.tx.isEmpty()) {
        QVector<double> y;
        for (const auto& symbol : CurrentOfdm.tx) y.append(symbol.real());
        QVector<double> scaledX;
        double scaleFactor = 1e6;  // 10^6

        scaledX.reserve(CurrentOfdm.t.size());
        for (double x : CurrentOfdm.t) {
            scaledX.append(x * scaleFactor);
        }
        auto xRange = plotsOfTimeDomain[0]->xAxis->range();
        plotsOfTimeDomain[0]->graph(0)->setData(scaledX, y);
        if(first_in) {
            plotsOfTimeDomain[0]->graph(0)->rescaleAxes();
            first_in = false;
        }
        else plotsOfTimeDomain[0]->xAxis->setRange(xRange);
        plotsOfTimeDomain[0]->replot();
    }
}

void GraphPlotting::PlotPSDPlots(const QVector<QVector<double>>& PSDs, const QVector<QVector<double>>& freqs)
{
    if(!PSDs.isEmpty() && !freqs.isEmpty()) {
        plotsOfPSD[0]->graph(0)->setData(freqs[0], PSDs[0]);
        plotsOfPSD[1]->graph(0)->setData(freqs[0], PSDs[1]);
        plotsOfPSD[0]->graph(0)->rescaleAxes();
        plotsOfPSD[1]->graph(0)->rescaleAxes();
        plotsOfPSD[0]->replot();
        plotsOfPSD[1]->replot();

        plotsOfPSD[3]->graph(0)->setData(freqs[0], PSDs[0]);
        plotsOfPSD[3]->graph(1)->setData(freqs[0], PSDs[1]);
        plotsOfPSD[3]->graph(0)->rescaleAxes();
        plotsOfPSD[3]->replot();
    }
}

void GraphPlotting::PlotPaCurve(PaCurve& PACurve, const QList<QAction*> actions)
{
    QString xLabel, yLabel;
    x.resize(PACurve.point_num);
    y.resize(PACurve.point_num);

    if (actions[0]->isChecked() && actions[3]->isChecked()) {
        x = PACurve.P_in_norm.linear;
        y     = PACurve.P_out_norm.linear;
        wp_x  = PACurve.Working_point_linear_norm.x;
        wp_y  = PACurve.Working_point_linear_norm.y;
        xLabel = "Pвх/Pнас, Вт";
        yLabel = "Pвх/Pнас, Вт";
    }
    else if (actions[0]->isChecked() && actions[2]->isChecked()) {
        x     = PACurve.P_in_abs.linear;
        y     = PACurve.P_out_abs.linear;
        wp_x  = PACurve.Working_point_linear_abs.x;
        wp_y  = PACurve.Working_point_linear_abs.y;
        xLabel = "Pвх, Вт";
        yLabel = "Pвых, Вт";
    }
    else if (actions[1]->isChecked() && actions[3]->isChecked()) {
        x     = PACurve.P_in_norm.dB;
        y     = PACurve.P_out_norm.dB;
        wp_x  = PACurve.Working_point_dB_norm.x;
        wp_y  = PACurve.Working_point_dB_norm.y;
        xLabel = "Pвх/Pнас, дБВт";
        yLabel = "Pвых/Pнас, дБВт";
    }
    else {
        x     = PACurve.P_in_abs.dB;
        y     = PACurve.P_out_abs.dB;
        wp_x  = PACurve.Working_point_dB_abs.x;
        wp_y  = PACurve.Working_point_dB_abs.y;
        xLabel = "Pвх, дБВт";
        yLabel = "Pвых, дБВт";
    }

    plotPaCurve->graph(0)->setData(x, y);
    plotPaCurve->graph(1)->setData(x, PACurve.Phi);
    plotPaCurve->graph(2)->setData(wp_x, wp_y);
    plotPaCurve->graph(3)->setData(wp_x, PACurve.Phi_work_grad);

    plotPaCurve->xAxis->setLabel(xLabel);
    plotPaCurve->yAxis->setLabel(yLabel);

    plotPaCurve->graph(0)->rescaleAxes();
    plotPaCurve->graph(1)->rescaleAxes();

    auto yRange = plotPaCurve->yAxis->range();
    double padding = (yRange.upper - yRange.lower) * 0.05;
    plotPaCurve->yAxis->setRange(yRange.lower, yRange.upper + padding);

    plotPaCurve->replot();
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
    PlotIdealSymConstellation("BPSK");
}

void GraphPlotting::InitializeTimeDomainPlotting(QVector<QWidget*> TimeDomainGraphWidgets)
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
        plotsOfTimeDomain[i]->replot();
        plotsOfTimeDomain[i]->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
        plotsOfTimeDomain[i]->axisRect()->setRangeZoom(Qt::Horizontal);
        plotsOfTimeDomain[i]->axisRect()->setRangeDrag(Qt::Horizontal);
    }
}

void GraphPlotting::InitializePSDPlotting(QVector<QWidget*> PSDGraphWidgets)
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
            plotsOfPSD[i]->graph(0)->setName("Перед PA");
            plotsOfPSD[i]->graph(1)->setName("После PA");
            plotsOfPSD[i]->legend->setVisible(true);
            plotsOfPSD[i]->graph(1)->setPen(QPen(Qt::red));
            plotsOfPSD[i]->graph(1)->setLineStyle(QCPGraph::lsLine);
            plotsOfPSD[i]->graph(1)->setScatterStyle(QCPScatterStyle::ssNone);
        }
        plotsOfPSD[i]->xAxis->setLabel("f, МГц");
        plotsOfPSD[i]->xAxis->setNumberFormat("f");
        plotsOfPSD[i]->xAxis->setNumberPrecision(1);
        plotsOfPSD[i]->yAxis->setLabel("Amplitude");
        plotsOfPSD[i]->replot();
        plotsOfPSD[i]->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
        plotsOfPSD[i]->axisRect()->setRangeZoom(Qt::Horizontal);
        plotsOfPSD[i]->axisRect()->setRangeDrag(Qt::Horizontal);
    }
}

void GraphPlotting::InitializeConstellationsPlotting(QVector<QWidget *> ConstellationsGraphWidgets)
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
        plotsOfConstellations[i]->replot();
    }
}
