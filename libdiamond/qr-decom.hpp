#pragma once
#ifndef DIAMOND_QR_DECOMPOSITION_HPP
#define DIAMOND_QR_DECOMPOSITION_HPP

#include "matrix.hpp"
#include "vector.hpp"
#include "util.hpp"

namespace Diamond {
namespace QR {

inline int Sgn(const double &x)
{
	return x < -EPS_DOUBLE ? -1 : 1;
}

inline int Sgn(const long double &x)
{
	return x < -EPS_LONG_DOUBLE ? -1 : 1;
}

/**
 * QR Decomposition without column pivoting
 * Input: matrix A
 * Return a std::pair of <Q, R> which satisfies A = Q * R
 */
template<typename _Td>
std::pair<Matrix<_Td>, Matrix<_Td>> QR_Decomposition(const Matrix<_Td> &a)
{
	Matrix<_Td> q(a.RowSize(), a.RowSize());
	Matrix<_Td> r(a);
	for (size_t i = 0; i < q.RowSize(); ++i) {
		q[i][i] = 1;
	}
	for (size_t k = 0; k < r.ColSize(); ++k) {
		_Td alpha = 0.0;
		for (size_t i = k; i < r.RowSize(); ++i) {
			alpha += r[i][k] * r[i][k];
		}
		alpha = -Sgn(r[k][k]) * sqrt(alpha);
		Vector<_Td> v(r.RowSize(), 0);
		for (size_t i = k; i < r.RowSize(); ++i) {
			v[i] = r[i][k];
		}
		v[k] -= alpha;
		_Td beta = 0;
		for (size_t i = 0; i < r.RowSize(); ++i) {
			beta += v[i] * v[i];
		}
		if (EqualZero(beta)) {
			continue;
		}
		for (size_t j = k; j < r.ColSize(); ++j) {
			_Td gama = 0;
			for (size_t i = 0; i < r.RowSize(); ++i) {
				gama += v[i] * r[i][j];
			}
			for (size_t i = 0; i < r.RowSize(); ++i) {
				r[i][j] = r[i][j] - (2 * gama) / beta * v[i];
			}
		}
		for (size_t j = 0; j < r.RowSize(); ++j) {
			_Td lambda = 0;
			for (size_t i = 0; i < r.RowSize(); ++i) {
				lambda += v[i] * q[j][i];
			}
			for (int i = 0; i < r.RowSize(); ++i) {
				q[j][i] = q[j][i] - (2 * lambda) / beta * v[i];
			}
		}
	}
	return std::make_pair(q, r);
}

/**
 * QR Decomposition with column pivoting
 * Input: a matrix A
 * Return <<Q, R>, P> which satisfies A = Q * R * P' where P' is the transpose of P.
 */
template<typename _Td>
std::pair<std::pair<Matrix<_Td>, Matrix<_Td>>, Matrix<_Td>> QR_DecompositionPivoting(const Matrix<_Td> &a)
{
	std::vector<size_t> P(a.ColSize());
	for (size_t i = 0; i < a.ColSize(); ++i) {
		P[i] = i;
	}
	Matrix<_Td> q(I<_Td>(a.ColSize())), r(a);
	for (size_t k = 0; k < r.ColSize(); ++k) {
		//----------Pivoting---------------
		size_t pivotCol = k;
		_Td pivotBeta = 0.0;
		for (size_t c = k; c < r.ColSize(); ++c) {
			_Td alpha = 0.0;
			for (size_t i = k; i < r.RowSize(); ++i) {
				alpha += r[i][c];
			}
			alpha = -Sgn(r[k][c]) * sqrt(alpha);
			Vector<_Td> v(r.RowSize(), 0);
			for (size_t i = k; i < r.RowSize(); ++i) {
				v[i] = r[i][c];
			}
			v[k] -= alpha;
			_Td beta = 0.0;
			for (size_t i = 0; i < r.RowSize(); ++i) {
				beta += v[i] * v[i];
			}
			if (fabs(beta) > fabs(pivotBeta)) {
				pivotBeta = beta;
				pivotCol = c;
			}
		}
		if (pivotCol != k) {
			for (size_t i = 0; i < r.RowSize(); ++i) {
				std::swap(r[i][k], r[i][pivotCol]);
			}
			std::swap(P[k], P[pivotCol]);
		}
		//----------Pivoting End-------------
		_Td alpha = 0.0;
		for (size_t i = k; i < r.RowSize(); ++i) {
			alpha += r[i][k] * r[i][k];
		}
		alpha = -Sgn(r[k][k]) * sqrt(alpha);
		Vector<_Td> v(r.RowSize(), 0);
		for (size_t i = k; i < r.RowSize(); ++i) {
			v[i] = r[i][k];
		}
		v[k] -= alpha;
		_Td beta = 0;
		for (size_t i = 0; i < r.RowSize(); ++i) {
			beta += v[i] * v[i];
		}
		if (EqualZero(beta)) {
			continue;
		}
		for (size_t j = k; j < r.ColSize(); ++j) {
			_Td gama = 0;
			for (size_t i = 0; i < r.RowSize(); ++i) {
				gama += v[i] * r[i][j];
			}
			for (size_t i = 0; i < r.RowSize(); ++i) {
				r[i][j] = r[i][j] - (2 * gama) / beta * v[i];
			}
		}
		for (size_t j = 0; j < r.RowSize(); ++j) {
			_Td lambda = 0;
			for (size_t i = 0; i < r.RowSize(); ++i) {
				lambda += v[i] * q[j][i];
			}
			for (int i = 0; i < r.RowSize(); ++i) {
				q[j][i] = q[j][i] - (2 * lambda) / beta * v[i];
			}
		}
	}
	Matrix<_Td> p(r.ColSize(), r.ColSize());
	for (size_t i = 0; i < r.ColSize(); ++i) {
		p[P[i]][i] = static_cast<_Td>(1);
	}
	return std::make_pair(std::make_pair(q, r), p);
}

}
}
#endif
