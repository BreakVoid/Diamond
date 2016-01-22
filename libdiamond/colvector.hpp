#pragma once

#include "vectorbase.hpp"

#ifndef DIAMOND_COL_VECTOR_HPP
#define DIAMOND_COL_VECTOR_HPP
namespace Diamond {
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
	Vector(const std::vector<_Td> &content)
		: VectorBase<_Td>::VectorBase(content)
	{}
	Vector(std::vector<_Td> &&content)
		: VectorBase<_Td>::VectorBase(content)
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
	Vector<_Td> & operator=(const Vector<_Td> &rhs)
	{
		this->data = rhs.data;
		this->size = rhs.size;
		return *this;
	}
	Vector<_Td> & operator=(Vector<_Td> &&rhs)
	{
		this->data = rhs.data;
		this->size = rhs.size;
		return *this;
	}
	friend std::ostream & operator<<(std::ostream &stream, const Vector<_Td> &vec)
	{
		std::ostream::fmtflags oldFlags = stream.flags();
		stream.precision(8);
		stream.setf(std::ios::fixed | std::ios::right);

		stream << '\n';
		for (size_t i = 0; i < vec.Size(); ++i) {
			stream << std::setw(15) << vec[i] << '\n';
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
Vector<_Td> operator*(const _Td &lhs, const Vector<_Td> &rhs)
{
	Vector<_Td> result(rhs.Size());
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
}

#endif