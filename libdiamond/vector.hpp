#pragma once

#ifndef DIAMOND_VECTOR_HPP
#define DIAMOND_VECTOR_HPP
#include "matrix.hpp"

namespace Diamond {
/**
 * Definition for vector base
 */
template<typename _Td>
class VectorBase {
protected:
	size_t size;
	std::vector<_Td> data;
public:
	size_t Size() const
	{
		return this->size;
	}
	const _Td & operator[](const size_t &pos) const
	{
		return this->data[pos];
	}
	_Td & operator[](const size_t &pos)
	{
		return this->data[pos];
	}
	VectorBase() : size(0) {}
	VectorBase(const size_t &_size, const _Td &initValue = static_cast<_Td>(0))
		: size(_size), data(size, initValue)
	{}
	VectorBase(const VectorBase<_Td> &vec)
		: size(vec.size), data(vec.data)
	{}
	VectorBase(VectorBase<_Td> &&vec)
		: size(vec.size), data(vec.data)
	{}
	VectorBase<_Td> &operator=(const VectorBase<_Td> &rhs)
	{
		this->size = rhs.size;
		this->data = rhs.data;
		return *this;
	}
	VectorBase<_Td> &operator=(VectorBase<_Td> &&rhs)
	{
		this->size = rhs.size;
		this->data = rhs.data;
		return *this;
	}
	virtual ~VectorBase() = default;
};

/**
 * Definition for column vector
 */
template<typename _Td>
class Vector : public VectorBase<_Td> {
public:
	Vector() : VectorBase<_Td>::VectorBase() {}
	Vector(const size_t &size, const _Td &initValue = static_cast<_Td>(0))
		: VectorBase<_Td>::VectorBase(size, initValue)
	{}
	Vector(const VectorBase<_Td> &vec)
		: VectorBase<_Td>::VectorBase(vec)
	{}
	Vector(VectorBase<_Td> &&vec)
		: VectorBase<_Td>::VectorBase(vec)
	{}
	Vector(const Vector<_Td> &vec)
		: VectorBase<_Td>::VectorBase(vec)
	{}
	Vector(Vector<_Td> &&vec)
		: VectorBase<_Td>::VectorBase(vec)
	{}
	Vector(const Matrix<_Td> &mat)
		: VectorBase<_Td>::VectorBase(mat.RowSize())
	{
		if (mat.ColSize() != 1) {
			throw std::invalid_argument("the matrix cannot be converted to a vector.");
		}
		for (size_t i = 0; i < this->size; ++i) {
			this->data[i] = mat[i][0];
		}
	}
	Vector(Matrix<_Td> &&mat)
		: VectorBase<_Td>::VectorBase(mat.RowSize())
	{
		if (mat.ColSize() != 1) {
			throw std::invalid_argument("the matrix cannot be converted to a vector.");
		}
		for (size_t i = 0; i < this->size; ++i) {
			this->data[i] = mat[i][0];
		}
	}
	friend std::ostream & operator<<(std::ostream &stream, const Vector<_Td> &vec)
	{
		std::ostream::fmtflags oldFlags = stream.flags();
		stream.precision(8);
		stream.setf(std::ios::fixed | std::ios::right);

		stream << '\n';
		for (size_t i = 0; i < vec.Size(); ++i) {
			stream << setw(15) << vec[i] << '\n';
		}
		stream.flags(oldFlags);
		return stream;
	}
};

template<typename _Td>
Vector<_Td> operator+(const Vector<_Td> &lhs, const Vector<_Td> &rhs)
{
	if (lhs.Size() != rhs.Size()) {
		throw std::invalid_argument("the sizes of two vectors cannot match.");
	}
	Vector<_Td> result(lhs.Size());
	for (size_t i = 0; i < result.Size(); ++i) {
		result[i] = lhs[i] + rhs[i];
	}
	return result;
}

template<typename _Td>
Vector<_Td> operator-(const Vector<_Td> &lhs, const Vector<_Td> &rhs)
{
	if (lhs.Size() != rhs.Size()) {
		throw std::invalid_argument("the sizes of two vectors cannot match.");
	}
	Vector<_Td> result(lhs.Size());
	for (size_t i = 0; i < result.Size(); ++i) {
		result[i] = lhs[i] - rhs[i];
	}
	return result;
}

template<typename _Td>
Vector<_Td> operator-(const Vector<_Td> &vec)
{
	Vector<_Td> result(vec.Size());
	for (size_t i = 0; i < result.Size(); ++i) {
		result[i] = -vec[i];
	}
	return result;
}

template<typename _Td>
Vector<_Td> operator-(Vector<_Td> &&vec)
{
	for (size_t i = 0; i < vec.Size(); ++i) {
		vec[i] = -vec[i];
	}
	return vec;
}

template<typename _Td>
Vector<_Td> operator*(const Vector<_Td> &lhs, const _Td &rhs)
{
	Vector<_Td> result(lhs.Size());
	for (size_t i = 0; i < result.Size(); ++i) {
		result[i] = lhs[i] * rhs;
	}
	return result;
}

template<typename _Td>
Vector<_Td> operator*(const _Td &hs, const Vector<_Td> &rhs)
{
	Vector<_Td> result(lhs.Size());
	for (size_t i = 0; i < result.Size(); ++i) {
		result[i] = lhs * rhs[i];
	}
	return result;
}

template<typename _Td>
Vector<_Td> operator/(const Vector<_Td> &lhs, const _Td &rhs)
{
	Vector<_Td> result(lhs.Size());
	for (size_t i = 0; i < result.Size(); ++i) {
		result[i] = lhs[i] / rhs;
	}
	return result;
}

/**
 * Definition for row vector
 */
template<typename _Td>
class VectorT : public VectorBase<_Td> {
public:
	VectorT() : VectorBase<_Td>::VectorBase() {}
	VectorT(const size_t &size, const _Td &initValue = static_cast<_Td>(0))
		: VectorBase<_Td>::VectorBase(size, initValue)
	{}
	VectorT(const VectorBase<_Td> &vec)
		: VectorBase<_Td>::VectorBase(vec)
	{}
	VectorT(VectorBase<_Td> &&vec)
		: VectorBase<_Td>::VectorBase(vec)
	{}
	VectorT(const VectorT<_Td> &vec)
		: VectorBase<_Td>::VectorBase(vec)
	{}
	VectorT(VectorT<_Td> &&vec)
		: VectorBase<_Td>::VectorBase(vec)
	{}
	VectorT(const Matrix<_Td> &mat)
		: VectorBase<_Td>::VectorBase(mat.ColSize())
	{
		if (mat.RowSize() != 1) {
			throw std::invalid_argument("the matrix cannot be converted to a vector.");
		}
		for (size_t i = 0; i < this->size; ++i) {
			this->data[i] = mat[0][i];
		}
	}
	VectorT(Matrix<_Td> &&mat)
		: VectorBase<_Td>::VectorBase(mat.RowSize())
	{
		if (mat.ColSize() != 1) {
			throw std::invalid_argument("the matrix cannot be converted to a vector.");
		}
		for (size_t i = 0; i < this->size; ++i) {
			this->data[i] = mat[i][0];
		}
	}
	friend std::ostream & operator<<(std::ostream &stream, const VectorT<_Td> &vec)
	{
		std::ostream::fmtflags oldFlags = stream.flags();
		stream.precision(8);
		stream.setf(std::ios::fixed | std::ios::right);

		stream << '\n';
		for (size_t i = 0; i < vec.Size(); ++i) {
			stream << setw(15) << vec[i];
		}
		stream << '\n';
		stream.flags(oldFlags);
		return stream;
	}
};

template<typename _Td>
VectorT<_Td> operator+(const VectorT<_Td> &lhs, const VectorT<_Td> &rhs)
{
	if (lhs.Size() != rhs.Size()) {
		throw std::invalid_argument("the sizes of two vectors cannot match.");
	}
	VectorT<_Td> result(lhs.Size());
	for (size_t i = 0; i < result.Size(); ++i) {
		result[i] = lhs[i] + rhs[i];
	}
	return result;
}

template<typename _Td>
VectorT<_Td> operator-(const VectorT<_Td> &lhs, const VectorT<_Td> &rhs)
{
	if (lhs.Size() != rhs.Size()) {
		throw std::invalid_argument("the sizes of two vectors cannot match.");
	}
	VectorT<_Td> result(lhs.Size());
	for (size_t i = 0; i < result.Size(); ++i) {
		result[i] = lhs[i] - rhs[i];
	}
	return result;
}

template<typename _Td>
VectorT<_Td> operator-(const VectorT<_Td> &vec)
{
	VectorT<_Td> result(vec.Size());
	for (size_t i = 0; i < result.Size(); ++i) {
		result[i] = -vec[i];
	}
	return result;
}

template<typename _Td>
VectorT<_Td> operator-(VectorT<_Td> &&vec)
{
	for (size_t i = 0; i < vec.Size(); ++i) {
		vec[i] = -vec[i];
	}
	return vec;
}

template<typename _Td>
VectorT<_Td> operator*(const VectorT<_Td> &lhs, const _Td &rhs)
{
	VectorT<_Td> result(lhs.Size());
	for (size_t i = 0; i < result.Size(); ++i) {
		result[i] = lhs[i] * rhs;
	}
	return result;
}

template<typename _Td>
VectorT<_Td> operator*(const _Td &hs, const VectorT<_Td> &rhs)
{
	VectorT<_Td> result(lhs.Size());
	for (size_t i = 0; i < result.Size(); ++i) {
		result[i] = lhs * rhs[i];
	}
	return result;
}

template<typename _Td>
VectorT<_Td> operator/(const VectorT<_Td> &lhs, const _Td &rhs)
{
	VectorT<_Td> result(lhs.Size());
	for (size_t i = 0; i < result.Size(); ++i) {
		result[i] = lhs[i] / rhs;
	}
	return result;
}

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
}

#endif

