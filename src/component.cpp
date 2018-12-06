// component.cpp
//
// Copyright (c) 2018 Ben Lindsay <benjlindsay@gmail.com>

#include "component.hpp"

void Component::init_site_grid_vars() {
  int n_weights = sim->n_subgrid_points;
  // int(std::pow(sim->mesh_order + 1, sim->dim));
  site_grid_indices = ArrayXXiR(n_sites, n_weights);
  site_grid_weights = ArrayXXdR(n_sites, n_weights);
}

void Component::calculate_site_grid_weights(int i_site,
                                            ArrayXi &subgrid_center_indices) {
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
  // grid_utils::get_spline_weights(dx_from_nearest_grid_point, grid_weights,
  // sim);
  double *grid_weights_i =
      site_grid_weights.data() + i_site * sim->dim * (sim->mesh_order + 1);
  sim->get_spline_weights(dx_from_nearest_grid_point, grid_weights_i);
}

void Component::add_site_to_grid(int i_site, ArrayXi &subgrid_center_indices) {
  int left_shift = sim->mesh_order / 2 + sim->mesh_order % 2;
  ArrayXXi subgrid_indices = sim->weight_subgrid_index_shifts.rowwise() +
                             subgrid_center_indices.transpose() - left_shift;
  int *subgrid_indices_ptr = subgrid_indices.data();
  int *Nx_ptr = sim->Nx.data();
  int *xyz_shifts_ptr = new int[sim->dim];
  for (int i = 0; i < sim->n_subgrid_points; i++) {
    int ix_global, iy_global, iz_global, global_index;
    ix_global = subgrid_indices_ptr[i];
    ix_global = (ix_global + Nx_ptr[0]) % Nx_ptr[0];
    iy_global = subgrid_indices_ptr[i + sim->n_subgrid_points];
    iy_global = (iy_global + Nx_ptr[1]) % Nx_ptr[1];
    if (sim->dim == 2) {
      global_index = sim->get_global_index(ix_global, iy_global);
    } else {
      iz_global = subgrid_indices_ptr[i + 2 * sim->n_subgrid_points];
      iz_global = (iz_global + Nx_ptr[2]) % Nx_ptr[2];
      global_index = sim->get_global_index(ix_global, iy_global, iz_global);
    }
    if (global_index < 0 || global_index >= sim->M) {
      utils::die("Bad index!");
    }

    int *weight_subgrid_index_shifts_ptr =
        sim->weight_subgrid_index_shifts.data();
    for (int d = 0; d < sim->dim; d++) {
      xyz_shifts_ptr[d] =
          weight_subgrid_index_shifts_ptr[d * sim->n_subgrid_points + i];
    }
    // Calculate the weight to add to the grid
    double total_weight = 1.0;
    double *grid_weights_i_site_ptr =
        site_grid_weights.data() + i_site * sim->dim * (sim->mesh_order + 1);
    for (int d = 0; d < sim->dim; d++) {
      total_weight *= grid_weights_i_site_ptr[d * (sim->mesh_order + 1) +
                                              xyz_shifts_ptr[d]];
    }
    total_weight /= sim->grid_point_volume;
    // Actually add site to grid
    Species_Type species = static_cast<Species_Type>(site_types[i_site]);
    rho_center_map[species][global_index] += total_weight;
  }
  delete[] xyz_shifts_ptr;
}

void Component::calculate_grid_densities() {
  // Zero out rho_center fields
  for (auto it = rho_center_map.begin(); it != rho_center_map.end(); it++) {
    // it->second points to the actual rho_center density array, for example
    // rho_center_map[A]
    it->second = ArrayXd::Zero(sim->ML);
  }
  ArrayXi subgrid_center_indices(sim->dim);
  for (int i_site = 0; i_site < n_sites; i_site++) {
    // Fill grid_weights array so that grid_weights(i, j) contains the weight
    // for dimension i of parameter j, where 0 <= j <= mesh_order + 1
    calculate_site_grid_weights(i_site, subgrid_center_indices);
    add_site_to_grid(i_site, subgrid_center_indices);
  }
}

void Component::move_particles() {
  ArrayXd diffusion_coeffs(n_sites);
  for (int i = 0; i < n_sites; i++) {
    Species_Type species = static_cast<Species_Type>(site_types[i]);
    diffusion_coeffs[i] = sim->diffusion_coeff_map[species];
  }
  site_coords += site_forces.colwise() * diffusion_coeffs * sim->timestep;
  utils::enforce_pbc(site_coords, sim->Lx);
}