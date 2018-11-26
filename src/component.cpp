// component.cpp
//
// Copyright (c) 2018 Ben Lindsay <benjlindsay@gmail.com>

#include "component.hpp"

void Component::calculate_site_grid_weights(int i_site,
                                            ArrayXi &subgrid_center_indices,
                                            ArrayXXd &grid_weights) {
  int mesh_order = sim->mesh_order;
  int dim = sim->dim;
  ArrayXd dx = sim->dx;
  ArrayXi Nx = sim->Nx;
  ArrayXd Lx = sim->Lx;

  ArrayXd site_coords_i = site_coords.row(i_site);
  ArrayXd dx_from_nearest_grid_point(dim);
  if (mesh_order % 2 == 0) {
    // Calculate distance to nearest grid points if order is even
    subgrid_center_indices = ((site_coords_i + 0.5 * dx) / dx).cast<int>();
    // If coordinates are near the positive wall in any dimension, make
    // the
    // nearest grid point 0 (periodic boundary conditions)
    for (int d = 0; d < dim; d++) {
      if (subgrid_center_indices[d] >= Nx[d]) {
        subgrid_center_indices[d] -= Nx[d];
        site_coords_i[d] -= Lx[d];
      }
    }
    dx_from_nearest_grid_point =
        site_coords_i - dx * subgrid_center_indices.cast<double>();
  } else {
    // Calculate distance to nearest mid-point between grid points if order
    // is odd
    subgrid_center_indices = (site_coords_i / dx).cast<int>();
    dx_from_nearest_grid_point =
        site_coords_i - dx * (subgrid_center_indices.cast<double>() + 0.5);
    // TODO: Ask about corresponsing code from original, why check for >=
    // Nx[j]? why use boundary grid point instead of 0 grid point?
  }
  grid_utils::get_spline_weights(dx_from_nearest_grid_point, grid_weights, sim);
  if (sim->iter == 0 && i_site == 0) {
    std::cout << grid_weights << std::endl;
  }
}

void Component::add_site_to_grid(int i_site, ArrayXi &subgrid_center_indices,
                                 ArrayXXd &grid_weights) {
  int left_shift = sim->mesh_order / 2 + sim->mesh_order % 2;
  ArrayXXi subgrid_indices =
      subgrid_center_indices.transpose().replicate(sim->n_subgrid_points, 1) -
      left_shift + sim->weight_subgrid_index_shifts;
  for (int i = 0; i < sim->n_subgrid_points; i++) {
    int ix_global, iy_global, iz_global, global_index;
    ix_global = subgrid_indices(i, 0);
    ix_global = (ix_global + sim->Nx[0]) % sim->Nx[0];
    iy_global = subgrid_indices(i, 1);
    iy_global = (iy_global + sim->Nx[1]) % sim->Nx[1];
    if (sim->dim == 2) {
      global_index = sim->get_global_index(ix_global, iy_global);
    } else {
      iz_global = subgrid_indices(i, 2);
      iz_global = (iz_global + sim->Nx[2]) % sim->Nx[2];
      global_index = sim->get_global_index(ix_global, iy_global, iz_global);
    }
    if (global_index < 0 || global_index >= sim->M) {
      utils::die("Bad index!");
    }
    int x_shift, y_shift, z_shift;
    x_shift = sim->weight_subgrid_index_shifts(i, 0);
    y_shift = sim->weight_subgrid_index_shifts(i, 1);
    // Calculate the weight to add to the grid
    double total_weight;
    if (sim->dim == 2) {
      total_weight = grid_weights(0, x_shift) * grid_weights(1, y_shift) /
                     sim->grid_point_volume;
    } else {
      z_shift = sim->weight_subgrid_index_shifts(i, 2);
      total_weight = grid_weights(0, x_shift) * grid_weights(1, y_shift) *
                     grid_weights(2, z_shift) / sim->grid_point_volume;
    }
    double x_weight = grid_weights(0, x_shift);
    double y_weight = grid_weights(1, y_shift);
    if (x_weight < 0) {
      utils::die("x_weight < 0");
    }
    if (y_weight < 0) {
      utils::die("y_weight < 0");
    }
    // Actually add site to grid
    Species_Type species = static_cast<Species_Type>(site_types[i_site]);
    rho_center_map[species][global_index] += total_weight;
  }
}

void Component::calculate_grid_densities() {
  ArrayXi subgrid_center_indices(sim->dim);
  ArrayXXd grid_weights(sim->dim, sim->mesh_order + 1);
  for (int i_site = 0; i_site < n_sites; i_site++) {
    // Fill grid_weights array so that grid_weights(i, j) contains the weight
    // for dimension i of parameter j, where 0 <= j <= mesh_order + 1
    calculate_site_grid_weights(i_site, subgrid_center_indices, grid_weights);
    add_site_to_grid(i_site, subgrid_center_indices, grid_weights);
  }
}