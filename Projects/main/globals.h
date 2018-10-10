#pragma once

#include <Partio.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCore>

#include <iostream>
#include <cmath>
#include <stdexcept>

constexpr int dim = 3;
using T = double;

#define CELL_SIZE 1.f

using vec2i = Eigen::Matrix<int, 2, 1>;
using vec3i = Eigen::Matrix<int, 3, 1>;

using vec3 = Eigen::Matrix<T, 3, 1>;
using mat3 = Eigen::Matrix<T, 3, 3>;

using vecXD = Eigen::VectorXd;

const T PI = 3.14159265358979323846264338327950288;

double floorX(double a) {
    return std::floor(a);
}

vec3 floor(vec3 v) {
	vec3 ret = v.unaryExpr(&floorX);
	return ret;
}

mat3 floor(mat3 m) {
	mat3 ret = m.unaryExpr(&floorX);
	return ret;
}

vec3 sign(vec3 v) {
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

T rng() { return ((T) std::rand() / (RAND_MAX)); }