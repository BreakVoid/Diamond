#pragma once
#ifndef DIAMOND_INTERPOLATION_HPP
#define DIAMOND_INTERPOLATION_HPP

#include "matrix.hpp"
#include "vector.hpp"
#include "util.hpp"
#include "lu-decom.hpp"

#include <vector>
#include <map>
#include <stdexcept>

namespace Diamond {
namespace Interpolation {

template<typename _Td>
inline _Td Sqr(const _Td &x)
{
	return x * x;
}

template<typename _Td>
inline _Td Cub(const _Td &x)
{
	return x * x * x;
}

template<typename _Td>
std::vector<Vector<_Td>> shootingMethodOnce(const std::vector<std::pair<_Td, _Td>> &data, const _Td &startSndDeri)
{
	if (data.size() <= 1) {
		throw std::invalid_argument("Need more points when interpolating.");
	}
	std::vector<Vector<_Td>> result;
	/**
	 * Initial segment, the first and second rows are used to 
	 * make sure the spline function pass though the data point,
	 * and the third and fourth rows are used to guarantee f' and f''
	 * equal to the given value.
	 */
	Matrix<_Td> firstA(4, 4);
	Vector<_Td> firstB(4);
	//--- first row ---
	firstA[0][0] = static_cast<_Td>(1);
	firstA[0][1] = data[0].first;
	firstA[0][2] = Sqr(data[0].first);
	firstA[0][3] = Cub(data[0].first);
	firstB[0] = data[0].second;
	//--- second row ---
	firstA[1][0] = static_cast<_Td>(1);
	firstA[1][1] = data[1].first;
	firstA[1][2] = Sqr(data[1].first);
	firstA[1][3] = Cub(data[1].first);
	firstB[1] = data[1].second;
	//--- third row ---
	firstA[2][0] = static_cast<_Td>(0);
	firstA[2][1] = static_cast<_Td>(1);
	firstA[2][2] = static_cast<_Td>(2) * data[0].first;
	firstA[2][3] = static_cast<_Td>(3) * Sqr(data[0].first);
	firstB[2] = startSndDeri;
	//--- fourth row ---
	firstA[3][0] = firstA[3][1] = static_cast<_Td>(0);
	firstA[3][2] = static_cast<_Td>(2);
	firstA[3][3] = static_cast<_Td>(6) * data[0].first;
	firstB[3] = 0;
	result.push_back(LU::SolveLinearEquationSystem(firstA, firstB));
	for (size_t i = 1; i < data.size() - 1; ++i) {
		/* The same to the matrix above. */
		Matrix<_Td> A(4, 4);
		Vector<_Td> b(4);
		//--- first row ---
		A[0][0] = static_cast<_Td>(1);
		A[0][1] = data[i].first;
		A[0][2] = Sqr(data[i].first);
		A[0][3] = Cub(data[i].first);
		b[0] = data[i].second;
		//--- second row ---
		A[1][0] = static_cast<_Td>(1);
		A[1][1] = data[i + 1].first;
		A[1][2] = Sqr(data[i + 1].first);
		A[1][3] = Cub(data[i + 1].first);
		b[1] = data[i + 1].second;
		//--- third row ---
		A[2][0] = static_cast<_Td>(0);
		A[2][1] = static_cast<_Td>(1);
		A[2][2] = static_cast<_Td>(2) * data[i].first;
		A[2][3] = static_cast<_Td>(3) * Sqr(data[i].first);
		b[2] = result.back()[1] + 2 * result.back()[2] * data[i].first + 3 * result.back()[3] * Sqr(data[i].first);
		//--- fourth row ---
		A[3][0] = A[3][1] = static_cast<_Td>(0);
		A[3][2] = static_cast<_Td>(2);
		A[3][3] = static_cast<_Td>(6) * data[i].first;
		b[3] = static_cast<_Td>(2) * result.back()[2] + static_cast<_Td>(6) * result.back()[3] * data[i].first;
		result.push_back(LU::SolveLinearEquationSystem(A, b));
	}
	return result;
}

template<typename _Td>
std::vector<std::vector<_Td>> NaturalSplineInterpolationEco(const std::vector<std::pair<_Td, _Td>> &data)
{
	if (data.size() <= 1) {
		throw std::invalid_argument("Need more data points.");
	}
	Vector<_Td> r = GenerateRandomVector<_Td>(2);
	std::vector<std::vector<Vector<_Td>>> res = {
		shootingMethodOnce(data, r[0]),
		shootingMethodOnce(data, r[1])
	};
	_Td dd1 = 2 * res[0].back()[2] + 6 * res[0].back()[3] * data.back().first;
	_Td dd2 = 2 * res[1].back()[2] + 6 * res[1].back()[3] * data.back().first;
	_Td lambda = -dd2 / (dd1 - dd2);
	std::vector<std::vector<_Td>> finalRes;
	for (size_t i = 0; i < data.size() - 1; ++i) {
		Vector<_Td> cur = res[0][i] * lambda + res[1][i] * (1 - lambda);
		std::vector<_Td> v(4);
		for (size_t j = 0; j < cur.Size(); ++j) {
			v[j] = cur[j];
		}
		finalRes.push_back(v);
	}
	return finalRes;
}

}
}

#endif /* end of DIAMOND_INTERPOLATION_HPP */