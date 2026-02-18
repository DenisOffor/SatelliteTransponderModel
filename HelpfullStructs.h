#ifndef HELPFULLSTRUCTS_H
#define HELPFULLSTRUCTS_H
#include <QVector>
#include <qmath.h>
#include <QPair>
#include <QString>
#include <QRandomGenerator>
#include <complex>
#include <cmath>
#include <random>


struct Curve2D {
    QVector<double> dB;
    QVector<double> linear;

    void resize(int n) {
        dB.resize(n);
        linear.resize(n);
    }
};

struct Point {
    QVector<double> x;
    QVector<double> y;

    void resize(int n) {
        x.resize(n);
        y.resize(n);
    }
};

struct Symbols{
    QVector<std::complex<double>> tr_sym_clean;
    QVector<std::complex<double>> tr_sym_noisy;
    QVector<std::complex<double>> rec_sym_noisy;
    QVector<int> data;

    void resize(int n) {
        data.resize(n);
        tr_sym_clean.resize(n);
        tr_sym_noisy.resize(n);
        rec_sym_noisy.resize(n);
    }
};

struct PaCurve {
    int point_num = 0;

    Curve2D P_in_abs;
    Curve2D P_in_norm;
    Curve2D P_out_abs;
    Curve2D P_out_norm;

    Point Working_point_dB_abs;
    Point Working_point_linear_abs;
    Point Working_point_dB_norm;
    Point Working_point_linear_norm;

    QVector<double> Phi;
    QVector<double> Phi_work_grad;
    QVector<double> r_in;
    QVector<double> Aout;

    explicit PaCurve(int size = 1000) {
        resize(size);
    }

    void resize(int size) {
        point_num = size;
        P_in_abs.resize(size);
        P_in_norm.resize(size);
        P_out_abs.resize(size);
        P_out_norm.resize(size);

        Phi.resize(size);
        r_in.resize(size);
        Aout.resize(size);

        Working_point_dB_abs.resize(1);
        Working_point_linear_abs.resize(1);
        Working_point_dB_norm.resize(1);
        Working_point_linear_norm.resize(1);

        Phi_work_grad.resize(1);
    }
};

struct Source {
    QString ModType;
    int NumSym;
    int SNRSymdB;
    int M;

    QString SigType;

    int SC_f_carrier;
    int SC_symrate;
    double SC_rolloff;
    int SC_filter_length;
    QString SC_FilterType;

    int OFDM_f_carrier;
    int OFDM_Nfft;
    int OFDM_GB_DC;
    int OFDM_GB_Nyq;
    int OFDM_cycle_prefix;

    int FDMA_f_carrier;
    int FDMA_symrate;
    int FDMA_num_subcarriers;
    int FDMA_step_carrier;

    int oversampling;
    int fs;
    int SNRSig;
};

struct NeedToRecalc{
    bool RecalcSymbols;
    bool RecalcNoiseSym;
    bool RecalcNoiseSig;
    bool RecalcOfdmFc;
    bool RecalcOfdmSig;
    bool FullRecalc;

    void init() {
        RecalcSymbols = true;
        RecalcNoiseSym = true;
        RecalcNoiseSig = true;
        RecalcOfdmFc = true;
        RecalcOfdmSig = true;
        FullRecalc = true;
    }

    void clear() {
        RecalcSymbols = false;
        RecalcNoiseSym = false;
        RecalcNoiseSig = false;
        RecalcOfdmFc = false;
        RecalcOfdmSig = false;
        FullRecalc = false;
    }
};

#endif // HELPFULLSTRUCTS_H
