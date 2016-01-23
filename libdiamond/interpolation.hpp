#pragma once
#ifndef DIAMOND_INTERPOLATION_HPP
#define DIAMOND_INTERPOLATION_HPP

#include "matrix.hpp"
#include "vector.hpp"
#include "util.hpp"
#include "lu-decom.hpp"
#include "polynomial.hpp"

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

/**
 * use shooting method to give a natural spline interpolation for the given data points.
 * The spline funcitons used are all polynomials of the degree 3.
 * If the input data has n data points, it will return n - 1 polynomials with 4-tuple.
 */
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

/**
 * a class to get Bspline interpolation.
 */
template<typename _Td>
class Bspline {
private:
	std::vector<std::pair<_Td, _Td>> data;
	Vector<_Td> coff;
public:
	size_t Index(const long long &i)
	{
		return 5 + i;
	}
	inline _Td T(const size_t &index)
	{
		return data[index].first;
	}
	inline _Td Y(const size_t &index)
	{
		return data[index].second;
	}
private:
	inline _Td v(const long long &k, const long long &i, const _Td &t)
	{
		return (t - T(Index(i))) / (T(Index(i + k) - T(Index(i))));
	}
	inline Polynomial<_Td> v(const long long &k, const long long &i)
	{
		return Polynomial<_Td>({
			static_cast<_Td>(1) / (T(Index(i + k)) - T(Index(i))),
			-T(Index(i)) / (T(Index(i + k)) - T(Index(i)))
		});
	}
	_Td B(const long long &k, const long long &i, const _Td &t)
	{
		if (k == 0) {
			if (T(Index(i)) <= t && t < T(Index(i + 1))) {
				return static_cast<_Td>(1);
			} else {
				return static_cast<_Td>(0);
			}
		} else {
			return v(k, i, t) * B(k - 1, i, t) + (1 - v(k, i + 1, t)) * B(k - 1, i + 1, t);
		}
	}
	Polynomial<_Td> B_Poly(const long long &k, const long long &i, const _Td &t)
	{
		if (k == 0) {
			if (T(Index(i)) <= t && t <= T(Index(i + 1))) {
				return Polynomial<_Td>(1);
			} else {
				return Polynomial<_Td>(0);
			}
		} else {
			return v(k, i) * B_Poly(k - 1, i, t) + (1 - v(k, i + 1)) * B_Poly(k - 1, i + 1, t);
		}
	}
public:
	Bspline(const std::vector<std::pair<_Td, _Td>> &origin) : data()
	{
		if (origin.size() <= 1) {
			throw std::invalid_argument("Need more data points to make an interpolation.");
		}
		const auto n = origin.size();
		for (_Td i = static_cast<_Td>(-5); i < static_cast<_Td>(0); i += static_cast<_Td>(1)) {
			data.push_back(std::make_pair((origin[1].first - origin[0].first) * i + origin[0].first, static_cast<_Td>(0)));
		}
		for (size_t i = 0; i < origin.size(); ++i) {
			data.push_back(origin[i]);
		}
		for (_Td i = static_cast<_Td>(1); i < static_cast<_Td>(4); i += static_cast<_Td>(1)) {
			data.push_back(std::make_pair((origin[n - 1].first - origin[n - 2].first) * i + origin[n - 1].first, static_cast<_Td>(0)));
		}
		Matrix<_Td> A(n + 2, n + 2);
		Vector<_Td> b(n + 2);
		for (size_t i = 0; i < n; ++i) {
			A[i][i] = B(3, i - 3, T(Index(i)));
			A[i][i + 1] = B(3, i - 2, T(Index(i)));
			A[i][i + 2] = B(3, i - 1, T(Index(i)));
			std::cerr << A[i][i] << '\t' << A[i][i + 1] << '\t' << A[i][i + 2] << std::endl;
			b[i] = Y(Index(i));
		}
		A[n][0] = Derivate(B_Poly(3, -3, T(Index(0))), 2)(T(Index(0)));
		A[n][1] = Derivate(B_Poly(3, -2, T(Index(0))), 2)(T(Index(0)));
		A[n][2] = Derivate(B_Poly(3, -1, T(Index(0))), 2)(T(Index(0)));
		std::cerr << A[n][0] << '\t' << A[n][1] << '\t' << A[n][2] << '\t' << std::endl;
		b[n] = 0;
		A[n + 1][n - 1] = Derivate(B_Poly(3, n - 4, T(Index(n - 1))), 2)(T(Index(n - 1)));
		A[n + 1][n] = Derivate(B_Poly(3, n - 3, T(Index(n - 1))), 2)(T(Index(n - 1)));
		A[n + 1][n + 1] = Derivate(B_Poly(3, n - 2, T(Index(n - 1))), 2)(T(Index(n - 1)));
		std::cerr << A[n + 1][n - 1] << '\t' << A[n + 1][n] << '\t' << A[n + 1][n + 1] << std::endl;
		b[n + 1] = 0;
		coff = LU::SolveLinearEquationSystem(A, b);
	}
	_Td y(const _Td &t)
	{
		_Td res = static_cast<_Td>(0);
		for (size_t i = 0; i < coff.Size(); ++i) {
			res += coff[i] * B(3, i - 3, t);
		}
		return res;
	}
};

};
}

#endif /* end of DIAMOND_INTERPOLATION_HPP */