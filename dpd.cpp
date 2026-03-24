#include "dpd.h"
#include <cmath>

VectorXd DPD::DPDsolve_least_squares(const MatrixXd& A, const VectorXd& b) {
    // Нахождение решения МНК через нормальное уравнение
    return (A.transpose() * A).ldlt().solve(A.transpose() * b);
}
