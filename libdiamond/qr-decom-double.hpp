#pragma once

#ifndef DIAMOND_QR_HPP
#define DIAMOND_QR_HPP

#include "matrix.hpp"
#include "vector.hpp"

namespace Diamond {
namespace QR_double {

extern const double EPS;
int sgn(const double &x);
std::pair<Matrix<double>, Matrix<double>> QR_Decomposition(const Matrix<double> &a);
std::pair<std::pair<Matrix<double>, Matrix<double>>, Matrix<double>> QR_DecompositionPivoting(const Matrix<double> &a);

}
}

#endif