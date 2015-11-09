#pragma once

#ifndef DIAMOND_VECTOR_HPP
#define DIAMOND_VECTOR_HPP
#include "matrix.hpp"

namespace Diamond {

template<typename _Td>
class Vector : public Matrix<_Td> {
public:
	Vector() : Matrix<_Td>::Matrix() {}
	Vector(const size_t &_n_rows)
		: Matrix<_Td>::Matrix(_n_rows, 1) {}
	Vector(const size_t &_n_rows, const _Td &fillValue)
		: Matrix<_Td>::Matrix(_n_rows, 1, fillValue) {}
	inline const size_t & Size() const
	{
		return this->RowSize();
	}
	_Td & operator[](const size_t &pos)
	{
		return this->data[pos][0];
	}
	const _Td & operator[](const size_t &pos) const
	{
		return this->data[pos][0];
	}
};

template<typename _Td>
class VectorT : public Matrix<_Td> {
public:
	VectorT() : Matrix<_Td>::Matrix() {}
	VectorT(const size_t &_n_cols)
		: Matrix<_Td>::Matrix(1, _n_cols) {}
	VectorT(const size_t &_n_cols, const _Td &fillValue)
		: Matrix<_Td>::Matrix(1, _n_cols, fillValue) {}
	inline const size_t & Size() const
	{
		return this->ColSize();
	}
	_Td & operator[](const size_t pos)
	{
		return this->data[0][pos];
	}
	const _Td & operator[](const size_t &pos) const
	{
		return this->data[0][pos];
	}
};

template<typename _Td>
Vector<_Td> Transpose(const VectorT<_Td> &vT)
{
	Vector<_Td> res(vT.Size());
	for (size_t i = 0; i < vT.Size(); ++i) {
		res[i] = vT[i];
	}
	return res;
}

template<typename _Td>
VectorT<_Td> Transpose(const Vector<_Td> &v)
{
	VectorT<_Td> res(v.Size());
	for (size_t i = 0; i < v.Size(); ++i) {
		res[i] = v[i];
	}
	return res;
}

template<typename _Td>
Vector<_Td> operator*(const _Td &b, const Vector<_Td> v)
{
	Vector<_Td> res(v);
	for (size_t i = 0; i < v.Size(); ++i) {
		res[i] = b * v[i];
	}
	return res;
}

template<typename _Td>
Vector<_Td> operator*(const Vector<_Td> v, const _Td &b)
{
	Vector<_Td> res(v);
	for (size_t i = 0; i < v.Size(); ++i) {
		res[i] = b * v[i];
	}
	return res;
}

template<typename _Td>
Vector<_Td> operator/(const Vector<_Td> v, const _Td &b)
{
	Vector<_Td> res(v);
	for (size_t i = 0; i < v.Size(); ++i) {
		res[i] = v[i] / b;
	}
	return res;
}

template<typename _Td>
Vector<_Td> operator+(const Vector<_Td> v1, const Vector<_Td> v2)
{
	if (v1.Size() != v2.Size()) {
		//throw std::exception("Cannot add two vectors of different sizes.");
	}
	Vector<_Td> res(v1.Size());
	for (size_t i = 0; i < v1.Size(); ++i) {
		res[i] = v1[i] + v2[i];
	}
	return res;
}

template<typename _Td>
VectorT<_Td> operator*(const _Td &b, const VectorT<_Td> v)
{
	VectorT<_Td> res(v);
	for (size_t i = 0; i < v.Size(); ++i) {
		res[i] = b * v[i];
	}
	return res;
}

template<typename _Td>
VectorT<_Td> operator*(const VectorT<_Td> v, const _Td &b)
{
	VectorT<_Td> res(v);
	for (size_t i = 0; i < v.Size(); ++i) {
		res[i] = b * v[i];
	}
	return res;
}

template<typename _Td>
VectorT<_Td> operator/(const VectorT<_Td> v, const _Td &b)
{
	VectorT<_Td> res(v);
	for (size_t i = 0; i < v.Size(); ++i) {
		res[i] = v[i] / b;
	}
	return res;
}

template<typename _Td>
VectorT<_Td> operator+(const VectorT<_Td> v1, const VectorT<_Td> v2)
{
	if (v1.Size() != v2.Size()) {
		throw std::exception("Cannot add two vectors of different sizes.");
	}
	VectorT<_Td> res(v1.Size());
	for (size_t i = 0; i < v1.Size(); ++i) {
		res[i] = v1[i] + v2[i];
	}
	return res;
}


template<typename _Td>
Matrix<_Td> operator*(const Vector<_Td> &v, const VectorT<_Td> &vT)
{
	Matrix<_Td> c(v.Size(), vT.Size(), 0);
	for (size_t i = 0; i < v.Size(); ++i) {
		for (size_t j = 0; j < vT.Size(); ++j) {
			c[i][j] = v[i] * vT[j];
		}
	}
	return c;
}

template<typename _Td>
_Td operator*(const VectorT<_Td> &vT, const Vector<_Td> &v)
{
	if (vT.Size() != v.Size()) {
		throw std::exception("Cannot calculate the interal product of vectors");
	}
	_Td c = 0;
	for (size_t i = 0; i < vT.Size(); ++i) {
		c += vT[i] * v[i];
	}
	return c;
}

}

#endif

