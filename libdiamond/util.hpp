#pragma once
#ifndef DIAMOND_UTIL_HPP
#define DIAMOND_UTIL_HPP

#include "matrix.hpp"
#include "vector.hpp"
#include <complex>
#include <random>
#include <ctime>

namespace Diamond {

static const double EPS_DOUBLE = 1e-8;
static const long double EPS_LONG_DOUBLE = 1e-10;
static std::default_random_engine engine(time(NULL));

inline bool EqualZero(const double &x)
{
	return -EPS_DOUBLE < x && x < EPS_DOUBLE;
}

inline bool EqualZero(const long double &x)
{
	return -EPS_LONG_DOUBLE < x && x < EPS_LONG_DOUBLE;
}

template<typename _Td>
bool EqualZero(const std::complex<_Td> &x)
{
	return EqualZero(x.real()) && EqualZero(x.imag());
}

template<typename _Td>
Vector<_Td> GenerateRandomVector(const size_t &n, const _Td &minValue = static_cast<_Td>(0), const _Td &maxValue = static_cast<_Td>(1))
{
	Vector<_Td> res(n);
	std::uniform_real_distribution<_Td> random(minValue, maxValue);
	for (size_t i = 0; i < n; ++i) {
		res[i] = random(engine);
	}
	return res;
}

template<typename _Td>
Matrix<_Td> GenerateRandomMatrix(const size_t &n, const size_t &m, const _Td &minValue = static_cast<_Td>(0), const _Td &maxValue = static_cast<_Td>(1))
{
	Matrix<_Td> res(n, m);
	std::uniform_real_distribution<_Td> random(minValue, maxValue);
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < m; ++j) {
			res[i][j] = random(engine);
		}
	}
	return res;
}

template<typename _Td>
Matrix<_Td> GenerateRandomSymmetricMatrix(const size_t &n, const _Td &minValue = static_cast<_Td>(0), const _Td &maxValue = static_cast<_Td>(1))
{
	Matrix<_Td> res(n, n);
	std::uniform_real_distribution<_Td> random(minValue, maxValue);
	for (size_t i = 0; i < n; ++i) {
		res[i][i] = random(engine);
		for (size_t j = i + 1; j < n; ++j) {
			res[i][j] = res[j][i] = random(engine);
		}
	}
	return res;
}


}

#endif