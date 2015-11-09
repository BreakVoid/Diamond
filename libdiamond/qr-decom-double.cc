#include "qr-decom-double.hpp"
#include <map>
#include <cmath>
#include <iostream>

namespace Diamond {
namespace QR_double {

extern const double EPS = 1e-9;

int sgn(const double &x)
{
	if (x < -EPS) {
		return -1;
	} else {
		return 1;
	}
}

std::pair<Matrix<double>, Matrix<double>> QR_Decomposition(const Matrix<double> &a)
{
	Matrix<double> q(a.RowSize(), a.RowSize());
	Matrix<double> r(a);
	for (size_t i = 0; i < q.RowSize(); ++i) {
		q[i][i] = 1;
	}
	for (size_t k = 0; k < r.ColSize(); ++k) {
		double alpha = 0.0;
		for (size_t i = k; i < r.RowSize(); ++i) {
			alpha += r[i][k] * r[i][k];
		}
		alpha = -sgn(r[k][k]) * sqrt(alpha);
		Vector<double> v(r.RowSize(), 0);
		for (size_t i = k; i < r.RowSize(); ++i) {
			v[i] = r[i][k];
		}
		v[k] -= alpha;
		double beta = 0;
		for (size_t i = 0; i < r.RowSize(); ++i) {
			beta += v[i] * v[i];
		}
		if (-EPS < beta && beta < EPS) {
			continue;
		}
		for (size_t j = k; j < r.ColSize(); ++j) {
			double gama = 0;
			for (size_t i = 0; i < r.RowSize(); ++i) {
				gama += v[i] * r[i][j];
			}
			for (size_t i = 0; i < r.RowSize(); ++i) {
				r[i][j] = r[i][j] - (2 * gama) / beta * v[i];
			}
		}
		for (size_t j = 0; j < r.RowSize(); ++j) {
			double lambda = 0;
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

std::pair<std::pair<Matrix<double>, Matrix<double>>, Matrix<double>> QR_DecompositionPivoting(const Matrix<double> &a)
{
	std::vector<size_t> P(a.ColSize());
	for (size_t i = 0; i < a.ColSize(); ++i) {
		P[i] = i;
	}
	Matrix<double> q(I<double>(a.ColSize())), r(a);
	for (size_t k = 0; k < r.ColSize(); ++k) {
		//----------Pivoting---------------
		size_t pivotCol = k;
		double pivotBeta = 0.0;
		for (size_t c = k; c < r.ColSize(); ++c) {
			double alpha = 0.0;
			for (size_t i = k; i < r.RowSize(); ++i) {
				alpha += r[i][c];
			}
			alpha = -sgn(r[k][c]) * sqrt(alpha);
			Vector<double> v(r.RowSize(), 0);
			for (size_t i = k; i < r.RowSize(); ++i) {
				v[i] = r[i][c];
			}
			v[k] -= alpha;
			double beta = 0.0;
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
		double alpha = 0.0;
		for (size_t i = k; i < r.RowSize(); ++i) {
			alpha += r[i][k] * r[i][k];
		}
		alpha = -sgn(r[k][k]) * sqrt(alpha);
		Vector<double> v(r.RowSize(), 0);
		for (size_t i = k; i < r.RowSize(); ++i) {
			v[i] = r[i][k];
		}
		v[k] -= alpha;
		double beta = 0;
		for (size_t i = 0; i < r.RowSize(); ++i) {
			beta += v[i] * v[i];
		}
		if (-EPS < beta && beta < EPS) {
			continue;
		}
		for (size_t j = k; j < r.ColSize(); ++j) {
			double gama = 0;
			for (size_t i = 0; i < r.RowSize(); ++i) {
				gama += v[i] * r[i][j];
			}
			for (size_t i = 0; i < r.RowSize(); ++i) {
				r[i][j] = r[i][j] - (2 * gama) / beta * v[i];
			}
		}
		for (size_t j = 0; j < r.RowSize(); ++j) {
			double lambda = 0;
			for (size_t i = 0; i < r.RowSize(); ++i) {
				lambda += v[i] * q[j][i];
			}
			for (int i = 0; i < r.RowSize(); ++i) {
				q[j][i] = q[j][i] - (2 * lambda) / beta * v[i];
			}
		}
	}
	Matrix<double> p(r.ColSize(), r.ColSize());
	for (size_t i = 0; i < r.ColSize(); ++i) {
		p[P[i]][i] = 1;
	}
	return std::make_pair(std::make_pair(q, r), p);
}

}
}