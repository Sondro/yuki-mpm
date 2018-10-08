#pragma once
#include "globals.h"

template <typename T>
class GridData {
public:
	// Constructors 
	GridData(vec3i dims);
	GridData(const GridData &other); // copy

	GridData(GridData &&other) : GridData(other.dims) {
		swap(*this, other);
	} // move

	virtual ~GridData() {} // destructor

	GridData &operator=(GridData other) {
		swap(*this, other);
	} // copy assignment operator

	virtual T operator()(T i, T j, T k) = 0;
	virtual T operator()(vec3i idx) {
		return (*this)(idx[0], idx[1], idx[2]);
	}

	virtual T interpolate(vec3 pos);

	friend void swap(GridData &o1, GridData &o2) {
		using std::swap;
		swap(o1.mData, o2.mData);
		swap(o1.dims, o2.dims);
	}

	int flat(int x, int y, int z) const {
		return y * dims[0] + x + z * dims[0] * dims[1];
	}

	std::vector<T> mData;
	vec3i dims;
};

template <typename T>
class FaceData : public GridData<T> {
	virtual T operator()(T i, T j, T k) {
		return GridData<T>::mData[flat(std::ceil(i), std::ceil(j), std::ceil(k))];
	}
public:
};

template <typename T>
class CenterData : public GridData<T> {
public:
};