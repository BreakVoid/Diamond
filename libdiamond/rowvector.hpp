#pragma once

#include "vectorbase.hpp"

#ifndef DIAMOND_ROW_VECTOR_HPP
#define DIAMOND_ROW_VECTOR_HPP

namespace Diamond {
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
	VectorT(const std::vector<_Td> &content)
		: VectorBase<_Td>::VectorBase(content)
	{}
	VectorT(std::vector<_Td> &&content)
		: VectorBase<_Td>::VectorBase(content)
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
	VectorT<_Td> & operator=(const VectorT<_Td> &rhs)
	{
		this->data = rhs.data;
		this->size = rhs.size;
		return *this;
	}
	VectorT<_Td> & operator=(VectorT<_Td> &&rhs)
	{
		this->data = rhs.data;
		this->size = rhs.size;
		return *this;
	}
	friend std::ostream & operator<<(std::ostream &stream, const VectorT<_Td> &vec)
	{
		std::ostream::fmtflags oldFlags = stream.flags();
		stream.precision(8);
		stream.setf(std::ios::fixed | std::ios::right);

		stream << '\n';
		for (size_t i = 0; i < vec.Size(); ++i) {
			stream << std::setw(15) << vec[i];
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
VectorT<_Td> operator*(const _Td &lhs, const VectorT<_Td> &rhs)
{
	VectorT<_Td> result(rhs.Size());
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
}
#endif