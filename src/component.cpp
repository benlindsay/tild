// component.cpp
//
// Copyright (c) 2018 Ben Lindsay <benjlindsay@gmail.com>

#include "component.hpp"

void Component::calculate_site_grid_weights(int i_site,
                                            ArrayXXd &grid_weights) {
  int mesh_order = sim->mesh_order;
  int dim = sim->dim;
  ArrayXd dx = sim->dx;
  ArrayXi Nx = sim->Nx;

  ArrayXd site_coords_i = site_coords.row(i_site);
  ArrayXi nearest_grid_indices(dim);
  ArrayXd dx_from_nearest_grid_point(dim);
  if (mesh_order % 2 == 0) {
    // Calculate distance to nearest grid points if order is even
    nearest_grid_indices = ((site_coords_i + 0.5 * dx) / dx).cast<int>();
    // If coordinates are near the positive wall in any dimension, make
    // the
    // nearest grid point 0 (periodic boundary conditions)
    for (int d = 0; d < dim; d++) {
      if (nearest_grid_indices[d] >= Nx[d]) {
        nearest_grid_indices[d] -= Nx[d];
      }
    }
    dx_from_nearest_grid_point =
        site_coords_i - dx * nearest_grid_indices.cast<double>();
  } else {
    // Calculate distance to nearest mid-point between grid points if order
    // is odd
    nearest_grid_indices = (site_coords_i / dx).cast<int>();
    dx_from_nearest_grid_point =
        site_coords_i - dx * (nearest_grid_indices.cast<double>() + 0.5);
    // TODO: Ask about corresponsing code from original, why check for >=
    // Nx[j]? why use boundary grid point instead of 0 grid point?
  }
  grid_utils::get_spline_weights(dx_from_nearest_grid_point, grid_weights, sim);
  if (sim->iter == 0 && i_site == 0) {
    std::cout << grid_weights << std::endl;
  }
}

void Component::calculate_grid_densities() {
  ArrayXXd grid_weights(sim->dim, sim->mesh_order + 1);
  for (int i_site = 0; i_site < n_sites; i_site++) {
    // Fill grid_weights array so that grid_weights(i, j) contains the weight
    // for dimension i of parameter j, where 0 <= j <= mesh_order + 1
    calculate_site_grid_weights(i_site, grid_weights);
  }
}