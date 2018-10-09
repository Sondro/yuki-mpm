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