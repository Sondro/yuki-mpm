#pragma once
#include "globals.h"

/**
 * GridData stores data at the grid center. Indexing is done as such:
 * \f$x_{i, j, k}\f$ such that \f$i, j, k = 1, \cdots, N\f$
 */
template <typename U>
class GridData {
public:
	// Constructors 
	GridData(vec3i dims) :
		length((dims[0] + 1) * 
			  (dims[1] + 1) * 
			  (dims[2] + 1)),
		mData(length, U(0)),
		dims(dims) {}
	GridData(const GridData<U> &other) :
		mData(other.mData),
		dims(other.dims) {} // copy

	GridData(GridData<U> &&other) : GridData(other.dims) {
		swap(*this, other);
	} // move

	virtual ~GridData() {} // destructor

	GridData<U> &operator=(GridData<U> other) {
		swap(*this, other);
	} // copy assignment operator

	/**
	 * Returns a reference to the value located at the grid cell
	 * \f$x_{i, j, k}\f$. To retrieve an interpolated value, use
	 * interpolate(vec3 pos) instead.
	 *
	 * Throws an std::out_of_range exception if the provided
	 * grid indices are out of range of this grid.
	 */
	virtual U &operator()(int i, int j, int k) {
		int idx = 0;
		try {
			idx = flat(i, j, k);
		} catch (std::out_of_range &) {
			throw std::out_of_range("GridData<U>::operator() : indices out of range");
		}
		return mData[idx];
	}

	/**
	 * Returns the value located at the grid cell
	 * \f$x_{i, j, k}\f$. To retrieve an interpolated value, use
	 * interpolate(vec3 pos) instead.
	 */
	virtual U operator()(int i, int j, int k) const {
		int idx = 0;
		try {
			idx = flat(i, j, k);
		} catch (std::out_of_range &) {
			throw std::out_of_range("GridData<U>::operator() : indices out of range");
		}
		return mData[idx];
	}

	virtual U &operator()(vec3i idx) {
		return (*this)(idx[0], idx[1], idx[2]);
	}
	virtual U operator()(vec3i idx) const {
		return (*this)(idx[0], idx[1], idx[2]);
	}

	friend void swap(GridData<U> &o1, GridData<U> &o2) {
		using std::swap;
		swap(o1.mData, o2.mData);
		swap(o1.dims, o2.dims);
	}

	int flat(int x, int y, int z) const {
		if (x < 0 || y < 0 || z < 0 ||
			x > dims[0] - 1 ||
			y > dims[1] - 1 ||
			z > dims[1] - 1) {
			throw std::out_of_range("GridData<U>::flat(int, int, int) : indices out of range");
		}
		return x + y * dims[0] + z * dims[0] * dims[1];
	}

	/**
	 * Convert world space coordinates pos to grid space coordinates
	 */
	vec3 worldToGrid(const vec3 &pos) {
		if (pos[0] < 0 || pos[1] < 0 || pos[2] < 0 ||
			pos[0] > dims[0] * CELL_SIZE ||
			pos[1] > dims[1] * CELL_SIZE ||
			pos[2] > dims[2] * CELL_SIZE) {
			throw std::out_of_range("GridData<U>worldToSelf(const vec3 &) : position out of range");
		}
		return pos / CELL_SIZE;
	}

	/**
	 * Return the cell that a particular world space position is in.
	 * This index will correspond to the lower x-y-z corner.
	 */
	vec3i cellOf(const vec3 &pos) {
		vec3 gridSpace = worldToGrid(pos);
		return vec3i(static_cast<int>(gridSpace[0]),
					 static_cast<int>(gridSpace[1]),
					 static_cast<int>(gridSpace[2]));
	}

	/**
	 * Return the cell space coordinates of this point where the origin is
	 * the min corner of the cell the point resides in
	 */
	vec3 worldToCell(const vec3 &pos) {
		vec3i cell = cellOf(pos);
		vec3 gridSpace = worldToGrid(pos);
		return gridSpace - cell;
	}

    /**
     * Resets the data of this gridData to 0.
     */
    void reset() {
        for (auto it = mData.begin(); it != mData.end(); ++it) {
            *it = U(0);
        }
    }

	int length;
	std::vector<U> mData;
	vec3i dims;
};
