#pragma once

#ifndef DIAMOND_VECTOR_HPP
#define DIAMOND_VECTOR_HPP
#include "vectorbase.hpp"
#include "rowvector.hpp"
#include "colvector.hpp"
#include <cmath>

namespace Diamond {

/**
* the transpose of vectors
*/

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
Vector<_Td> Transpose(VectorT<_Td> &&vT)
{
	return Vector<_Td>(VectorBase<_Td>(vT));
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
VectorT<_Td> Transpose(Vector<_Td> &&v)
{
	return VectorT<_Td>(VectorBase<_Td>(v));
}

/**
* Multiplication between vector and matrix
*/
template<typename _Td>
Vector<_Td> operator*(const Matrix<_Td> &A, const Vector<_Td> &v)
{
	if (A.ColSize() != v.Size()) {
		throw std::invalid_argument("The column size of matrix is different from the size of vector.");
	}
	Vector<_Td> result(A.RowSize(), 0);
	for (size_t i = 0; i < A.RowSize(); ++i) {
		for (size_t j = 0; j < A.ColSize(); ++j) {
			result[i] += A[i][j] * v[j];
		}
	}
	return result;
}

template<typename _Td>
VectorT<_Td> operator*(const VectorT<_Td> &vT, const Matrix<_Td> &A)
{
	if (vT.Size() != A.RowSize()) {
		throw std::invalid_argument("The row size of matrix is different from the size of the transpose vector.");
	}
	VectorT<_Td> result(A.ColSize(), 0);
	for (size_t i = 0; i < A.ColSize(); ++i) {
		for (size_t j = 0; j < A.RowSize(); ++j) {
			result[i] += vT[j] * A[j][i];
		}
	}
	return result;
}

/**
* Multiplication between column vector and row vector
*/
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
		throw std::invalid_argument("different vector size.");
	}
	_Td c = 0;
	for (size_t i = 0; i < vT.Size(); ++i) {
		c += vT[i] * v[i];
	}
	return c;
}

template<typename _Td>
_Td Norm(const VectorBase<_Td> &v, const int &order = 2)
{
	if (order <= 0) {
		throw std::invalid_argument("The order of norm should be greater than 0.");
	}
	_Td result = static_cast<_Td>(0);
	for (size_t i = 0; i < v.Size(); ++i) {
		result += pow(v[i], order);
	}
	return pow(result, static_cast<_Td>(1.0) / order);
}

}

#endif

