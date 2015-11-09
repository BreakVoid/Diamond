#pragma once

#ifndef DIAMOND_MATRIX_HPP
#define DIAMOND_MATRIX_HPP

#include <iostream>
#include <vector>

namespace Diamond {

template<typename _Td>
class Matrix {
protected:
	size_t n_rows = 0;
	size_t n_cols = 0;
	std::vector<std::vector<_Td>> data;
	class RowProxy {
		std::vector<_Td> &row;
	public:
		RowProxy(std::vector<_Td> & _row) : row(_row) {}
		_Td & operator[](const size_t &pos)
		{
			return row[pos];
		}
	};
	class ConstRowProxy {
		const std::vector<_Td> &row;
	public:
		ConstRowProxy(const std::vector<_Td> &_row) : row(_row) {}
		const _Td & operator[](const size_t &pos) const
		{
			return row[pos];
		}
	};
public:
	Matrix() {};
	Matrix(const size_t &_n_rows, const size_t &_n_cols)
		: n_rows(_n_rows), n_cols(_n_cols), data(std::vector<std::vector<_Td>>(n_rows, std::vector<_Td>(n_cols))) {}
	Matrix(const size_t &_n_rows, const size_t &_n_cols, const _Td &fillValue)
		: n_rows(_n_rows), n_cols(_n_cols), data(std::vector<std::vector<_Td>>(n_rows, std::vector<_Td>(n_cols, fillValue))) {}
	Matrix(const Matrix<_Td> &mat)
		: n_rows(mat.n_rows), n_cols(mat.n_cols), data(mat.data) {}
	Matrix(Matrix<_Td> &&mat) noexcept
		: n_rows(mat.n_rows), n_cols(mat.n_cols), data(mat.data) {}
	inline const size_t & RowSize() const
	{
		return n_rows;
	}
	inline const size_t & ColSize() const
	{
		return n_cols;
	}
	RowProxy operator[](const size_t &Kth)
	{
		return RowProxy(this->data[Kth]);
	}
	const ConstRowProxy operator[](const size_t &Kth) const
	{
		return ConstRowProxy(this->data[Kth]);
	}
	~Matrix() {};
};

/**
 * Sum of two matrics.
 */
template<typename _Td>
Matrix<_Td> operator+(const Matrix<_Td> &a, const Matrix<_Td> &b)
{
	if (a.RowSize() != b.RowSize() || a.ColSize() != b.ColSize()) {
		throw std::exception("Two matrics cannot be added together.");
	}
	Matrix<_Td> c(a.RowSize(), a.ColSize());
	for (size_t i = 0; i < a.RowSize(); ++i) {
		for (size_t j = 0; j < a.ColSize(); ++j) {
			c[i][j] = a[i][j] + b[i][j];
		}
	}
	return c;
}

/**
 * Multiplication of two matrics.
 */
template<typename _Td>
Matrix<_Td> operator*(const Matrix<_Td> &a, const Matrix<_Td> &b)
{
	if (a.ColSize() != b.RowSize()) {
		throw std::exception("Two matrics cannot be multiplied together.");
	}
	Matrix<_Td> c(a.RowSize(), b.ColSize(), 0);
	for (size_t i = 0; i < a.RowSize(); ++i) {
		for (size_t j = 0; j < b.ColSize(); ++j) {
			for (size_t k = 0; k < a.ColSize(); ++k) {
				c[i][j] += a[i][k] * b[k][j];
			}
		}
	}
	return c;
}

/**
 * Operations between a number and a matrix;
 */
template<typename _Td>
Matrix<_Td> operator*(const Matrix<_Td> &a, const _Td &b)
{
	Matrix<_Td> c(a.RowSize(), a.ColSize());
	for (size_t i = 0; i < a.RowSize(); ++i) {
		for (size_t j = 0; j < a.ColSize(); ++j) {
			c[i][j] = a[i][j] * b;
		}
	}
	return c;
}

template<typename _Td>
Matrix<_Td> operator*(const _Td &b, const Matrix<_Td> &a)
{
	Matrix<_Td> c(a.RowSize(), a.ColSize());
	for (size_t i = 0; i < a.RowSize(); ++i) {
		for (size_t j = 0; j < ColSize(); ++j) {
			c[i][j] = a[i][j] * b;
		}
	}
	return c;
}

template<typename _Td>
Matrix<_Td> operator/(const Matrix<_Td> &a, const double &b)
{
	Matrix<_Td> c(a.RowSize(), a.ColSize());
	for (size_t i = 0; i < a.RowSize(); ++i) {
		for (size_t j = 0; j < a.ColSize(); ++j) {
			c[i][j] = a[i][j] / b;
		}
	}
	return c;
}

template<typename _Td>
Matrix<_Td> Transpose(const Matrix<_Td> &a)
{
	Matrix<_Td> res(a.ColSize(), a.RowSize());
	for (size_t i = 0; i < a.ColSize(); ++i) {
		for (size_t j = 0; j < a.RowSize(); ++j) {
			res[i][j] = a[j][i];
		}
	}
	return res;
}

template<typename _Td>
std::ostream & operator<<(std::ostream &stream, const Matrix<_Td> &mat)
{
	stream << '\n';
	for (size_t i = 0; i < mat.RowSize(); ++i) {
		for (size_t j = 0; j < mat.ColSize(); ++j) {
			stream << '\t' << mat[i][j];
		}
		stream << '\n';
	}
	return stream;
}

template<typename _Td>
Matrix<_Td> I(const size_t &n)
{
	Matrix<_Td> res(n, n, 0);
	for (size_t i = 0; i < n; ++i) {
		res[i][i] = static_cast<_Td>(1);
	}
	return res;
}

}
#endif

