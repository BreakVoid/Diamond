#pragma once
#ifndef DIAMOND_LU_DECOM_HPP
#define DIAMOND_LU_DECOM_HPP

#include "matrix.hpp"
#include "vector.hpp"
#include "util.hpp"
#include <map>

namespace Diamond {
namespace LU {

/**
 * LU decomposition without column pivoting
 * Input: a square matrix A
 * Return: <L, U> which satisfies A = L * U
 */
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


/**
 * LU decomposition without column pivoting
 * Input: a square matrix A
 * Return: <<L, U>, P> which satisfies A = L * U * P' where P' is the transpose of P.
 */
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
			std::swap(p[k], p[maxColumn]);
			for (size_t i = 0; i < n; ++i) {
				std::swap(beta[i][k], beta[i][maxColumn]);
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

template<typename _Td>
Vector<_Td> SolveLowerTriangleSystem(const Matrix<_Td> &a, const Vector<_Td> &b)
{
	if (a.RowSize() != a.ColSize()) {
		throw std::invalid_argument("Cannot solve a no-square system.");
	}
	if (a.RowSize() != b.Size()) {
		throw std::invalid_argument("The size of matrix and vector is not matched.");
	}
	const auto n = a.RowSize();
	Vector<_Td> res(n);
	for (size_t i = 0; i < n; ++i) {
		_Td sum = static_cast<_Td>(0);
		for (size_t j = 0; j < i; ++j) {
			sum += a[i][j] * res[j];
		}
		res[i] = (b[i] - sum) / a[i][i];
	}
	return res;
}

template<typename _Td>
Vector<_Td> SolveUpperTriangleSystem(const Matrix<_Td> &a, const Vector<_Td> &b)
{
	if (a.RowSize() != a.ColSize()) {
		throw std::invalid_argument("Cannot solve a no-square system.");
	}
	if (a.RowSize() != b.Size()) {
		throw std::invalid_argument("The size of matrix and vector is not matched.");
	}
	const auto n = a.RowSize();
	Vector<_Td> res(n);
	for (size_t i = n; i > 0; --i) {
		_Td sum = 0;
		for (size_t j = n; j > i; --j) {
			sum += a[i - 1][j - 1] * res[j - 1];
		}
		res[i - 1] = (b[i - 1] - sum) / a[i - 1][i - 1];
	}
	return res;
}

template<typename _Td>
Vector<_Td> SolveLinearEquationSystem(const Matrix<_Td> &a, const Vector<_Td> &b)
{
	const auto luRes = LU_DecompositionPivoting(a);
	const auto y = SolveLowerTriangleSystem(luRes.first.first, b);
	const auto x = SolveUpperTriangleSystem(luRes.first.second, y);
	return luRes.second * x;
}

}
}


#endif