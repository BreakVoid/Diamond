#pragma once
#ifndef DIAMOND_UTIL_HPP
#define DIAMOND_UTIL_HPP

#include <complex>

namespace Diamond {

extern const double EPS_DOUBLE;
extern const long double EPS_LONG_DOUBLE;

bool EqualZero(const double &x);
bool EqualZero(const long double &x);

template<typename _Td>
bool EqualZero(const std::complex<_Td> &x)
{
	return EqualZero(x.real()) && EqualZero(x.imag());
}

}

#endif