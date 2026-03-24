#ifndef DPD_H
#define DPD_H

#include <vector>
#include <complex>
#include <Eigen/Dense>
#include <QDebug>

using namespace Eigen;

class DPD {
public:

private:
    VectorXd DPDsolve_least_squares(const MatrixXd& A, const VectorXd& b);
};

#endif // DPD_H
