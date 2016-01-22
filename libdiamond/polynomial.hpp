#pragma once
#ifndef DIAMOND_POLYNOMIAL_HPP
#define DIAMOND_POLYNOMIAL_HPP value

#include "util.hpp"

#include <vector>
#include <iostream>
#include <sstream>
#include <functional>
#include <algorithm>

/**
 * A class of polynomial to make addtion, substraction, multiplication between polynomials easier.
 * And provide an API to get derivative of a polynomial.
 */

namespace Diamond {

template<typename _Td>
class Polynomial {
protected:
	std::vector<_Td> c;
public:
	Polynomial() : c()
	{}
	Polynomial(const _Td &_c) : c(1, _c)
	{}
	Polynomial(const std::vector<_Td> &_c) : c(_c)
	{}
	Polynomial(std::vector<_Td> &&_c) : c(_c)
	{}
	Polynomial(const Polynomial<_Td> &polynomial) : c(polynomial.c)
	{}
	Polynomial(Polynomial<_Td> &&polynomial) : c(polynomial.c)
	{}
	Polynomial<_Td> Trimed() const
	{
		size_t res = c.size();
		while (res > 1 && EqualZero(c[res - 1])) {
			--res;
		}
		auto cc = c;
		cc.resize(res);
		return Polynomial<_Td>(std::move(cc));
	}
	void trim()
	{
		size_t res = c.size();
		while (res > 1 && EqualZero(c[res - 1])) {
			--res;
		}
		c.resize(res);	
	}
	size_t Deg()
	{
		size_t res = c.size();
		while (res > 1 && EqualZero(c[res - 1])) {
			--res;
		}
		c.resize(res);
		return res - 1;
	}
	size_t Deg() const
	{
		size_t res = c.size();
		while (res > 1 && EqualZero(c[res - 1])) {
			--res;
		}
		return res - 1;
	}
	_Td & operator[](const size_t &pos)
	{
		if (pos >= c.size()) {
			c.resize(pos + 1, static_cast<_Td>(0));
		}
		return c[pos];
	}
	const _Td operator[](const size_t &pos) const
	{
		if (pos >= c.size()) {
			return static_cast<const _Td>(0);
		}
		return c[pos];
	}
	_Td & operator()(const _Td &x) const
	{
		_Td res = static_cast<_Td>(0);
		_Td powX = static_cast<_Td>(1);
		for (size_t i = 0; i < c.size(); ++i) {
			res += powX * c[i];
			powX *= x;
		}
		return res;
	}
	friend std::ostream &operator<<(std::ostream &stream, const Polynomial<_Td> &polynomial)
	{
		size_t deg = polynomial.Deg();
		// std::cerr << deg << std::endl;
		for (size_t i = deg + 1; i > 0; --i) {
			if (EqualZero(polynomial[i - 1])) {
				continue;
			}
			std::ostringstream oss;
			oss << polynomial[i - 1];
			if (i > 1) {
				oss << 'x';
			}
			if (i > 2) {
				oss << "^" << i - 1;
			}
			std::string str = oss.str();
			if (i < deg + 1 && str[0] != '-') {
				stream << "+";
			}
			stream << str;
		}
		return stream;
	}
	friend Polynomial<_Td> operator+(const Polynomial<_Td> &lhs, const Polynomial<_Td> &rhs)
	{
		Polynomial<_Td> res;
		size_t maxDeg = std::max(lhs.Deg(), rhs.Deg());
		for (size_t i = 0; i < maxDeg + 1; ++i) {
			res[i] = lhs[i] + rhs[i];
		}
		return res.Trimed();
	}
	friend Polynomial<_Td> operator-(const Polynomial<_Td> &lhs, const Polynomial<_Td> &rhs)
	{
		Polynomial<_Td> res;
		size_t maxDeg = std::max(lhs.Deg(), rhs.Deg());
		for (size_t i = 0; i < maxDeg + 1; ++i) {
			res[i] = lhs[i] - rhs[i];
		}
		return res.Trimed();
	}
	friend Polynomial<_Td> operator*(const Polynomial<_Td> &lhs, const Polynomial<_Td> &rhs)
	{
		Polynomial<_Td> res;
		const auto lhsDeg = lhs.Deg();
		const auto rhsDeg = rhs.Deg();
		for (size_t i = 0; i < lhsDeg + 1; ++i) {
			for (size_t j = 0; j < rhsDeg + 1; ++j) {
				res[i + j] += lhs[i] * rhs[j];
			}
		}
		return res.Trimed();
	}
	friend Polynomial<_Td> operator/(const Polynomial<_Td> &lhs, const _Td &rhs)
	{
		Polynomial<_Td> res;
		const auto deg = lhs.Deg();
		for (size_t i = 0; i < deg + 1; ++i) {
			res[i] = lhs[i] / rhs;
		}
		return res;
	}
};

}

#endif /* end of DIAMOND_POLYNOMIAL_HPP */