#pragma once
#ifndef DIAMOND_OPTIMIZATION_HPP
#define DIAMOND_OPTIMIZATION_HPP

#include "util.hpp"
#include "vector.hpp"
#include <functional>

namespace Diamond {
namespace Optimization {

static const double TAU_DOUBLE = (sqrt(5.0) - 1) / 2;
static const long double TAU_LONG_DOUBLE = (sqrt(5.0L) - 1) / 2;

template<typename _Tp>
_Tp GoldenSectionSearch(std::function<double(const _Tp &)> f, const _Tp &origin, const _Tp &direction)
{
	const double &TAU = TAU_DOUBLE;
	double a = 0, b = 1e10;
	double x = a + (1 - TAU) * (b - a);
	double y = a + TAU * (b - a);
	double fx = f(origin + x * direction);
	double fy = f(origin + y * direction);
	while (b - a > EPS_LONG_DOUBLE) {
		if (fx > fy) {
			a = x;
			x = y;
			fx = fy;
			y = a + TAU * (b - a);
			fy = f(origin + y * direction);
		} else {
			b = y;
			y = x;
			fy = fx;
			x = a + (1 - TAU) * (b - a);
			fx = f(origin + x * direction);
		}
	}
	return origin + a * direction;
}

template<typename _Tp>
_Tp GoldenSectionSearchHD(std::function<long double(const _Tp &)> f, const _Tp &origin, const _Tp &direction)
{
	const long double &TAU = TAU_LONG_DOUBLE;
	long double a = 0, b = 1e10L;
	long double x = a + (1.0L - TAU) * (b - a);
	long double y = a + TAU * (b - a);
	long double fx = f(origin + x * direction);
	long double fy = f(origin + y * direction);
	while (b - a > EPS_LONG_DOUBLE) {
		if (fx > fy) {
			a = x;
			x = y;
			fx = fy;
			y = a + TAU * (b - a);
			fy = f(origin + y * direction);
		} else {
			b = y;
			y = x;
			fy = fx;
			x = a + (1 - TAU) * (b - a);
			fx = f(origin + x * direction);
		}
	}
	return origin + a * direction;
}

template<typename _Tp>
Vector<_Tp> RandomGuessSolution(const size_t &sizeDimension, const size_t &randomTime, std::function<_Tp(const Vector<_Tp> &)> Func)
{
	Vector<_Tp> bestVector = GenerateRandomVector(sizeDimension, static_cast<_Tp>(-1e5), static_cast<_Tp>(1e5));
	double bestValue = Func(bestVector);
	for (size_t i = 0; i < randomTime; ++i) {
		Vector<_Tp> newVector = GenerateRandomVector(sizeDimension, static_cast<_Tp>(-1e5), static_cast<_Tp>(1e5));
		double newValue = Func(newVector);
		if (newValue < bestValue) {
			bestVector = newVector;
			bestValue = newValue;
		}
	}
	return bestVector;
}

template<typename _Tp>
Vector<_Tp> ConjugateGradientMethod(
	const size_t &sizeDimension,
	std::function<_Tp(const Vector<_Tp> &)> Func, // a function f: _Tp^n --> double
	std::function<Vector<_Tp>(const Vector<_Tp> &)> Gradient // a function to calculate the gradient of f
	)
{
	Vector<_Tp> x = RandomGuessSolution<_Tp>(sizeDimension, 10 * sizeDimension, Func);
	Vector<_Tp> g = Gradient(x);
	Vector<_Tp> s = -g;
	while (true) {
		Vector<_Tp> newX = GoldenSectionSearchHD<Vector<_Tp>>(Func, x, s);
		if (EqualZero(Func(x) - Func(newX))) {
			return newX;
		}
		Vector<_Tp> newG = Gradient(newX);
		_Tp beta = (Transpose(newG) * newG) / (Transpose(g) * g);
		s = -newG + beta * s;
		g = newG;
		x = newX;
	}
	return x;
}

template<typename _Tp>
Vector<_Tp> ConjugateGradientMethod(
	const size_t &sizeDimension,
	std::function<long double(const Vector<_Tp> &)> Func, // a function f: _Tp^n --> double
	std::function<Vector<_Tp>(const Vector<_Tp> &)> Gradient, // a function to calculate the gradient of f
	const Vector<_Tp> &firstSolution
	)
{
	Vector<_Tp> x = firstSolution;
	Vector<_Tp> g = Gradient(x);
	Vector<_Tp> s = -g;
	while (true) {
		Vector<_Tp> newX = GoldenSectionSearchHD<Vector<_Tp>>(Func, x, s);
		if (EqualZero(Func(x) - Func(newX))) {
			return newX;
		}
		Vector<_Tp> newG = Gradient(newX);
		_Tp beta = (Transpose(newG) * newG) / (Transpose(g) * g);
		s = -newG + beta * s;
		g = newG;
		x = newX;
	}
	return x;
}

}
}

#endif