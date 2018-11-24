// grid_utils.hpp
//
// Copyright (c) 2018 Ben Lindsay <benjlindsay@gmail.com>

#ifndef GRID_UTILS_HPP
#define GRID_UTILS_HPP

#include "Eigen/Dense"
#include "sim.hpp"

namespace grid_utils {
void get_spline_weights(ArrayXd &dx_from_nearest_grid_point,
                        ArrayXXd &grid_weights, Sim *sim);
}

#endif  // GRID_UTILS_HPP
