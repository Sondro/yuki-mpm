#pragma once
#include "globals.h"

template <typename T>
class GridData {
public:
	// Constructors 
	GridData(vec3i dims);
	GridData(const GridData &other); // copy

	GridData(GridData<T> &&other) : GridData(other.dims) {
		swap(*this, other);
	} // move

	virtual ~GridData() {} // destructor

	GridData<T> &operator=(GridData<T> other) {
		swap(*this, other);
	} // copy assignment operator

	virtual T operator()(T i, T j, T k) {
		return this->mData[flat(i, j, k)];
	}
	virtual T operator()(vec3i idx) {
		return (*this)(idx[0], idx[1], idx[2]);
	}

	virtual T interpolate(vec3 pos);

	friend void swap(GridData<T> &o1, GridData<T> &o2) {
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
