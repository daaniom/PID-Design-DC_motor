#ifndef EXTENDED_KALMAN_FILTER_H
#define EXTENDED_KALMAN_FILTER_H

#include <BasicLinearAlgebra.h>
#include "mecotron.h"


void PredictionUpdate(const Matrix<1> &u, Matrix<3> &xhat, Matrix<3,3> &Phat);
void CorrectionUpdate(const Matrix<1> &y, Matrix<3> &xhat, Matrix<3,3> &Phat, Matrix<1> &nu, Matrix<1,1> &S);

#endif // EXTENDED_KALMAN_FILTER_H
