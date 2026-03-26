#ifndef HELPFULLSTRUCTS_H
#define HELPFULLSTRUCTS_H
#include <qmath.h>
#include <QPair>
#include <QString>
#include <QRandomGenerator>
#include <complex>
#include <cmath>
#include <random>
#include "QDebug"
#include "fft.h"

struct Curve2D {
    std::vector<double> dB;
    std::vector<double> linear;

    void resize(int n) {
        dB.resize(n);
        linear.resize(n);
    }
};

struct Point {
    std::vector<double> x;
    std::vector<double> y;

    void resize(int n) {
        x.resize(n);
        y.resize(n);
    }
};

struct Symbols{
    std::vector<std::complex<double>> tr_sym_clean;
    std::vector<std::complex<double>> tr_sym_noisy;
    std::vector<std::complex<double>> sym_after_pa;
    std::vector<std::complex<double>> rec_sym_noisy;
    std::vector<std::complex<double>> rec_sym_noisy_with_DPD;
    std::vector<int> data_tx;
    std::vector<int> data_rx;
    std::vector<int> data_rx_with_DPD;

    void resize(int n) {
        data_tx.resize(n);
        data_rx.resize(n);
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

    std::vector<double> Phi;
    std::vector<double> Phi_work_grad;
    std::vector<double> r_in;
    std::vector<double> Aout;

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

    QString PAModel;
    int linear_gain_dB;
    int IBO_dB;
    std::vector<double> SalehCoeffs;
    std::vector<double> RappCoeffs;
    std::vector<double> GhorbaniCoeffs;
    std::vector<double> FIRCoeffs;
    QString StaticNonlinModel;

    int MP_M;
    int MP_P;
    QString NormalizationType;
};

struct NeedToRecalc{
    bool RecalcSymbols;
    bool RecalcNoiseSym;
    bool RecalcNoiseSig;
    bool RecalcSig;
    bool FullRecalc;
    bool PARecalc;
    bool TimePlotsRescale;
    bool PaCurveReplot;
    bool DPDRecalc;

    void init() {
        RecalcSymbols = true;
        RecalcNoiseSym = true;
        RecalcNoiseSig = true;
        RecalcSig = true;
        FullRecalc = true;
        PARecalc = true;
        TimePlotsRescale = true;
        PaCurveReplot = true;
        DPDRecalc = true;
    }

    void clear() {
        RecalcSymbols = false;
        RecalcNoiseSym = false;
        RecalcNoiseSig = false;
        RecalcSig = false;
        FullRecalc = false;
        PARecalc = false;
        TimePlotsRescale = false;
        PaCurveReplot = false;
        DPDRecalc = false;
    }
};

struct GlobalResults {
    std::vector<std::complex<double>> tx_sig;
    std::vector<std::complex<double>> tx_plus_dpd_sig;
    std::vector<std::complex<double>> pa_sig;
    std::vector<std::complex<double>> pa_plus_dpd_sig;
    std::vector<std::complex<double>> rec_sig;
    std::vector<double> time;

    double BER_noDPD;
    double BER_withDPD;
    double EVM_noDPD;
    double EVM_withDPD;
    double ACLR_noDPD;
    double ACLR_withDPD;

    void resize(int size) {
        tx_sig.resize(size);
        tx_plus_dpd_sig.resize(size);
        pa_sig.resize(size);
        pa_plus_dpd_sig.resize(size);
        rec_sig.resize(size);
    }

    void clear() {
        tx_sig.clear();
        pa_sig.clear();
        pa_plus_dpd_sig.clear();
        rec_sig.clear();
    }
};

struct PSDdata {
    std::vector<std::vector<double>> f;
    std::vector<std::vector<double>> PSD;

    void resize(int n) {
        f.resize(n);
        PSD.resize(n);
    }

    void clear() {
        f.clear();
        PSD.clear();
    }
};

#endif // HELPFULLSTRUCTS_H
