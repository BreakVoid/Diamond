#pragma once
#ifndef DIAMOND_LU_DECOM_HPP
#define DIAMOND_LU_DECOM_HPP

#include "matrix.hpp"
#include "vector.hpp"
#include "util.hpp"
#include <map>

namespace Diamond {
namespace LU {

template<typename _Td>
std::pair<Matrix<_Td>, Matrix<_Td>> LU_Decomposition(const Matrix<_Td> &a)
{
	if (a.RowSize() != a.ColSize()) {
		throw std::invalid_argument("Cannot apply LU decomposition on a matrix which is not square.");
	}
	const size_t n = a.RowSize();
	Matrix<_Td> alpha(n, n);
	Matrix<_Td> beta = a;
	for (size_t k = 0; k < n; ++k) {
		if (EqualZero(beta[k][k])) {
			break;
		}
		for (size_t i = k; i < n; ++i) {
			alpha[i][k] = beta[i][k] / beta[k][k];
		}
		for (size_t j = k + 1; j < n; ++j) {
			for (size_t i = k + 1; i < n; ++i) {
				beta[i][j] = beta[i][j] - alpha[i][k] * beta[k][j];
			}
		}
	}
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < i; ++j) {
			beta[i][j] = static_cast<_Td>(0);
		}
	}
	return std::make_pair(alpha, beta);
}

template<typename _Td>
std::pair<std::pair<Matrix<_Td>, Matrix<_Td>>, Matrix<_Td>> LU_DecompositionPivoting(const Matrix<_Td> &a)
{
	if (a.RowSize() != a.ColSize()) {
		throw std::invalid_argument("Cannot apply LU decomposition on a matrix which is not square.");
	}
	const size_t n = a.RowSize();
	Matrix<_Td> alpha(n, n);
	Matrix<_Td> beta = a;
	std::vector<size_t> p(n);
	for (size_t i = 0; i < n; ++i) {
		p[i] = i;
	}
	for (size_t k = 0; k < n; ++k) {
		size_t maxColumn = k;
		for (size_t i = k + 1; i < n; ++i) {
			if (abs(beta[k][i]) > abs(beta[k][maxColumn])) {
				maxColumn = i;
			}
		}
		if (maxColumn != k) {
			swap(p[k], p[maxColumn]);
			for (size_t i = 0; i < n; ++i) {
				swap(beta[i][k], beta[i][maxColumn]);
			}
		}
		if (EqualZero(beta[k][k])) {
			break;
		}
		for (size_t i = k; i < n; ++i) {
			alpha[i][k] = beta[i][k] / beta[k][k];
		}
		for (size_t j = k + 1; j < n; ++j) {
			for (size_t i = k + 1; i < n; ++i) {
				beta[i][j] = beta[i][j] - alpha[i][k] * beta[k][j];
			}
		}
	}
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < i; ++j) {
			beta[i][j] = static_cast<_Td>(0);
		}
	}
	Matrix<_Td> P(n, n);
	for (size_t i = 0; i < n; ++i) {
		P[p[i]][i] = static_cast<_Td>(1);
	}
	return std::make_pair(std::make_pair(alpha, beta), P);
}

}
}


#endif