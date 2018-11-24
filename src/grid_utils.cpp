// grid_utils.cpp
//
// Copyright (c) 2018 Ben Lindsay <benjlindsay@gmail.com>

#include "grid_utils.hpp"

void grid_utils::get_spline_weights(ArrayXd &dx_from_nearest_grid_point,
                                    ArrayXXd &grid_weights, Sim *sim) {
  ArrayXd dx_norm = dx_from_nearest_grid_point / sim->dx;

  int mesh_order = sim->mesh_order;
  double scale = double(sim->M) / sim->V;
  if (mesh_order == 0) {
    // TODO: Is this true? Why the scale when grid weights sum to 1 along each
    // dimension in all other cases?
    grid_weights = 1.0 * scale;
  } else if (mesh_order == 1) {
    grid_weights.col(0) = (1.0 - 2.0 * dx_norm) / 2.0;
    grid_weights.col(1) = (1.0 + 2.0 * dx_norm) / 2.0;
  } else if (mesh_order == 2) {
    ArrayXd dx_norm_2 = dx_norm * dx_norm;

    grid_weights.col(0) = (1.0 - 4.0 * dx_norm + 4.0 * dx_norm_2) / 8.0;
    grid_weights.col(1) = (3.0 - 4.0 * dx_norm_2) / 4.0;
    grid_weights.col(2) = (1.0 + 4.0 * dx_norm + 4.0 * dx_norm_2) / 8.0;
  } else if (mesh_order == 3) {
    ArrayXd dx_norm_2 = dx_norm * dx_norm;
    ArrayXd dx_norm_3 = dx_norm * dx_norm_2;

    grid_weights.col(0) =
        (1.0 - 6.0 * dx_norm + 12.0 * dx_norm_2 - 8.0 * dx_norm_3) / 48.0;
    grid_weights.col(1) =
        (23.0 - 30.0 * dx_norm - 12.0 * dx_norm_2 + 24.0 * dx_norm_3) / 48.0;
    grid_weights.col(2) =
        (23.0 + 30.0 * dx_norm - 12.0 * dx_norm_2 - 24.0 * dx_norm_3) / 48.0;
    grid_weights.col(3) =
        (1.0 + 6.0 * dx_norm + 12.0 * dx_norm_2 + 8.0 * dx_norm_3) / 48.0;
  } else if (mesh_order == 4) {
    ArrayXd dx_norm_2 = dx_norm * dx_norm;
    ArrayXd dx_norm_3 = dx_norm * dx_norm_2;
    ArrayXd dx_norm_4 = dx_norm * dx_norm_3;

    grid_weights.col(0) = (1.0 - 8.0 * dx_norm + 24.0 * dx_norm_2 -
                           32.0 * dx_norm_3 + 16.0 * dx_norm_4) /
                          384.0;
    grid_weights.col(1) = (19.0 - 44.0 * dx_norm + 24.0 * dx_norm_2 +
                           16.0 * dx_norm_3 - 16.0 * dx_norm_4) /
                          96.0;
    grid_weights.col(2) =
        (115.0 - 120.0 * dx_norm_2 + 48.0 * dx_norm_4) / 192.0;
    grid_weights.col(3) = (19.0 + 44.0 * dx_norm + 24.0 * dx_norm_2 -
                           16.0 * dx_norm_3 - 16.0 * dx_norm_4) /
                          96.0;
    grid_weights.col(4) = (1.0 + 8.0 * dx_norm + 24.0 * dx_norm_2 +
                           32.0 * dx_norm_3 + 16.0 * dx_norm_4) /
                          384.0;
  } else if (mesh_order == 5) {
    ArrayXd dx_norm_2 = dx_norm * dx_norm;
    ArrayXd dx_norm_3 = dx_norm * dx_norm_2;
    ArrayXd dx_norm_4 = dx_norm * dx_norm_3;
    ArrayXd dx_norm_5 = dx_norm * dx_norm_4;

    grid_weights.col(0) =
        (1.0 - 10.0 * dx_norm + 40.0 * dx_norm_2 - 80.0 * dx_norm_3 +
         80.0 * dx_norm_4 - 32.0 * dx_norm_5) /
        3840.0;
    grid_weights.col(1) =
        (237.0 - 750.0 * dx_norm + 840.0 * dx_norm_2 - 240.0 * dx_norm_3 -
         240.0 * dx_norm_4 + 160.0 * dx_norm_5) /
        3840.0;
    grid_weights.col(2) =
        (841.0 - 770.0 * dx_norm - 440.0 * dx_norm_2 + 560.0 * dx_norm_3 +
         80.0 * dx_norm_4 - 160.0 * dx_norm_5) /
        1920.0;
    grid_weights.col(3) =
        (841.0 + 770.0 * dx_norm - 440.0 * dx_norm_2 - 560.0 * dx_norm_3 +
         80.0 * dx_norm_4 + 160.0 * dx_norm_5) /
        1920.0;
    grid_weights.col(4) =
        (237.0 + 750.0 * dx_norm + 840.0 * dx_norm_2 + 240.00 * dx_norm_3 -
         240.0 * dx_norm_4 - 160.0 * dx_norm_5) /
        3840.0;
    grid_weights.col(5) =
        (1.0 + 10.0 * dx_norm + 40.0 * dx_norm_2 + 80.0 * dx_norm_3 +
         80.0 * dx_norm_4 + 32.0 * dx_norm_5) /
        3840.0;
  } else {
    utils::die("get_spline_weights not set up for this interpolation order!\n");
  }
}