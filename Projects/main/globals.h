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

constexpr T DENSITY = 2.0;
constexpr T VOLUME = 1.0;

// The length, width, and height of each cell
// Each cell is a cube
constexpr T CELL_SIZE = 0.3;
constexpr T INV_CELL_SIZE = 1.0 / CELL_SIZE;

// The extent of grid (grid goes from {0, 0, 0} to {X, Y, Z}
// Along with the cell size, determines how many cells should
// be in the grid
constexpr T X_SIZE = 10;
constexpr T Y_SIZE = 10;
constexpr T Z_SIZE = 10;

constexpr int X_CELL_COUNT = static_cast<int>(X_SIZE / CELL_SIZE);
constexpr int Y_CELL_COUNT = static_cast<int>(Y_SIZE / CELL_SIZE);
constexpr int Z_CELL_COUNT = static_cast<int>(Z_SIZE / CELL_SIZE);

#define DEBUG 1
//#define APIC
#define STRESS 1
#define GRAVITY 1

constexpr T dt = 0.0001;
constexpr T k = 1000;
constexpr T nu = .3;
constexpr T V0 = 1e-3;

using vec2i = Eigen::Matrix<int, 2, 1>;
using vec3i = Eigen::Matrix<int, 3, 1>;

using vec3 = Eigen::Matrix<T, 3, 1>;
using mat3 = Eigen::Matrix<T, 3, 3>;
using mat4 = Eigen::Matrix<T, 4, 4>;

using vecXD = Eigen::VectorXd;

constexpr T PI = 3.14159265358979323846264338327950288;

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

static vec3 g = []() {
	vec3 g;
	g << 0, -9.8, 0;
	return g;
}();

constexpr T EPSILON = 100.0 * dt;
constexpr T WEIGHT_CONSTRAINT = 1.0;
static vec3 ZERO_VECTOR = []() {
	return vec3::Zero();
}();
