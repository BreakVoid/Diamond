#pragma once

#ifndef DIAMOND_VECTOR_BASE_HPP
#define DIAMOND_VECTOR_BASE_HPP
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
	VectorBase(const std::vector<_Td> &content)
		: size(content.size()), data(content)
	{}
	VectorBase(std::vector<_Td> &&content)
		: size(content.size()), data(content)
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

}

#endif

