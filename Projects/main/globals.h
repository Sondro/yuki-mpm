#pragma once

#include <Partio.h>
#include "voro++.hh"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCore>

#include <iostream>
#include <cmath>
#include <stdexcept>

constexpr int dim = 3;
using T = double;

#define CELL_SIZE (T(1.0))
#define INV_CELL_SIZE (T(1.f / CELL_SIZE))
#define DEBUG 1
#define dt (T(0.000001))
static T k = 100000;
static T nu = .3;
static T V0 = 1e-3;

using vec2i = Eigen::Matrix<int, 2, 1>;
using vec3i = Eigen::Matrix<int, 3, 1>;

using vec3 = Eigen::Matrix<T, 3, 1>;
using mat3 = Eigen::Matrix<T, 3, 3>;

using vecXD = Eigen::VectorXd;

const T PI = 3.14159265358979323846264338327950288;

inline double floorX(double a) {
    return std::floor(a);
}

inline vec3 floor(vec3 v) {
	vec3 ret = v.unaryExpr(&floorX);
	return ret;
}

inline mat3 floor(mat3 m) {
	mat3 ret = m.unaryExpr(&floorX);
	return ret;
}

inline vec3 sign(vec3 v) {
	vec3 ret;
	ret[0] = v[0] < 0 ? -1 : 1;
	ret[1] = v[1] < 0 ? -1 : 1;
	ret[2] = v[2] < 0 ? -1 : 1;

	ret[0] = v[0] == 0 ? 0 : ret[0];
	ret[1] = v[1] == 0 ? 0 : ret[1];
	ret[2] = v[2] == 0 ? 0 : ret[2];
	return ret;
}

using PointList = std::vector<vec3>;

inline T rng() { return ((T) std::rand() / (RAND_MAX)); }

static mat3 D = ((1.0 / 3.0) * CELL_SIZE * CELL_SIZE * mat3::Identity());
static mat3 D_INV = D.inverse();