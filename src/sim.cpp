// sim.cpp
//
// Copyright (c) 2017 Ben Lindsay <benjlindsay@gmail.com>

#include "sim.hpp"

Sim::Sim(YAML::Node input) {
  utils::print_one_line("Initializing Sim");
  // Initialize dim
  if (!input["dim"]) {
    utils::die("dim keyword not found in input file");
  } else {
    dim = input["dim"].as<int>();
  }

  iter = 0;

  // Initialize max_iter
  if (!input["max_iter"]) {
    max_iter = 100;
  } else {
    max_iter = input["max_iter"].as<int>();
  }

  // Initialize rho_0
  if (!input["rho_0"]) {
    utils::die("rho_0 not found in input file");
  } else {
    rho_0 = input["rho_0"].as<double>();
  }

  // Initialize box/grid variables
  init_box_vars(input);

  // Initialize components
  init_component_list(input);

  // Initialize forces
}

Sim::~Sim() {
  for (int i = 0; i < output_list.size(); i++) {
    delete output_list[i];
  }
  for (int i = 0; i < component_list.size(); i++) {
    delete component_list[i];
  }
}

void Sim::init_box_vars(YAML::Node input) {
  // Initialize Lx, Nx, and dx
  Lx = ArrayXd(dim);
  Nx = ArrayXi(dim);
  dx = ArrayXd(dim);

  // 2 of the 3 must be input, and we can calculate the 3rd from the other 2.
  // To figure out which 2 were input, we'll start with the set of those 3 and
  // delete the corresponding string from the set until we're left with just 1
  std::set<std::string> keys = {"lx", "nx", "dx"};
  for (YAML::const_iterator it = input.begin(); it != input.end(); ++it) {
    if (it->first.as<std::string>() == "lx") {
      if (input["lx"].size() != size_t(dim)) {
        utils::die("Length of lx must equal dim (" + std::to_string(dim) + ")");
      }
      for (std::size_t i = 0; i < size_t(dim); i++) {
        Lx[i] = input["lx"][i].as<double>();
      }
      keys.erase("lx");
    } else if (it->first.as<std::string>() == "nx") {
      if (input["nx"].size() != size_t(dim)) {
        utils::die("Length of nx must equal dim (" + std::to_string(dim) + ")");
      }
      for (std::size_t i = 0; i < size_t(dim); i++) {
        Nx[i] = input["nx"][i].as<int>();
      }
      keys.erase("nx");
    } else if (it->first.as<std::string>() == "dx") {
      if (input["dx"].size() != size_t(dim)) {
        utils::die("Length of dx must equal dim (" + std::to_string(dim) + ")");
      }
      for (std::size_t i = 0; i < size_t(dim); i++) {
        dx[i] = input["dx"][i].as<double>();
      }
      keys.erase("dx");
    }
  }
  if (keys.size() != 1) {
    utils::die("You must define exactly two of the following: Lx, Nx, dx");
  }

  // Calculate the 3rd variable (Lx, Nx, or dx)
  std::string remaining_var = *(keys.begin());
  if (remaining_var == "lx") {
    Lx = Nx.cast<double>() * dx;
  } else if (remaining_var == "nx") {
    // Calculate Nx rounded to nearest odd int (Odd values work better with
    // FFTW)
    // First, truncate Lx/dx
    Nx = (Lx / dx).cast<int>();
    // Next, add one to any even value in Nx.
    for (int i = 0; i < dim; i++) {
      if (Nx[i] % 2 == 0) {
        Nx[i] += 1;
      }
    }
    // Recalculate dx with the same Lx and the calculated Nx
    dx = Lx / Nx.cast<double>();
  } else if (remaining_var == "dx") {
    dx = Lx / Nx.cast<double>();
  }

  // Assign aggregate values
  V = Lx.prod();
  M = Nx.prod();
  grid_point_volume = dx.prod();  // Also = V / double(M)
  ML = M;

  local_grid_indices = ArrayXXi::Zero(ML, dim);
  int n_reps = 1;
  for (int d = 0; d < dim; d++) {
    for (int i = 0; i < ML; i++) {
      local_grid_indices(i, d) = (i / n_reps) % Nx[d];
    }
    n_reps *= Nx[d];
  }

  local_grid_coords = ArrayXXd::Zero(ML, dim);
  for (int d = 0; d < dim; d++) {
    local_grid_coords.col(d) = local_grid_indices.col(d).cast<double>() * dx[d];
  }

  if (!input["mesh_order"]) {
    mesh_order = 1;
  } else {
    mesh_order = input["mesh_order"].as<int>();
    if (mesh_order < 0 || mesh_order > 5) {
      utils::die("mesh_order must be between 0 and 5 inclusive!");
    }
  }

  n_subgrid_points = std::pow(mesh_order + 1, dim);
  weight_subgrid_index_shifts = ArrayXXi::Zero(n_subgrid_points, dim);
  n_reps = 1;
  for (int d = 0; d < dim; d++) {
    for (int i = 0; i < n_subgrid_points; i++) {
      weight_subgrid_index_shifts(i, d) = (i / n_reps) % (mesh_order + 1);
    }
    n_reps *= mesh_order + 1;
  }
}

void Sim::init_component_list(YAML::Node input) {
  utils::print_one_line("Running Sim::init_component_list()");
  if (component_list.size() > 0) {
    utils::die("component_list has already been filled!");
  }
  YAML::Node components = input["components"];
  for (std::size_t i = 0; i < components.size(); i++) {
    Component *comp =
        Component_Factory::New_Component(this, input, components[i]);
    component_list.push_back(comp);
  }
}

void Sim::write_outputs() {
  for (size_t i = 0; i < output_list.size(); i++) {
    if (output_list[i]->is_time_to_write()) {
      output_list[i]->write();
    }
  }
}

void Sim::run() {
  utils::print_one_line("Running " + description);
  // Write outputs for initial state
  std::cout << "mesh_order: " << mesh_order << std::endl;
  for (iter = 0; iter < max_iter; iter++) {
    calculate_grid_densities();
    calculate_forces();
    write_outputs();
    move_particles();
  }
  write_outputs();
}

void Sim::calculate_grid_densities() {
  for (int i_comp = 0; i_comp < component_list.size(); i_comp++) {
    Component *comp = component_list[i_comp];
    for (auto it = comp->rho_center_map.begin();
         it != comp->rho_center_map.end(); it++) {
      // it->second points to the actual rho_center density array, for example
      // rho_center_map[A]
      it->second = ArrayXd::Zero(ML);
    }
    comp->calculate_grid_densities();
  }
}

int Sim::get_global_index(int ix_global, int iy_global) {
  return ix_global + Nx[0] * iy_global;
}

int Sim::get_global_index(int ix_global, int iy_global, int iz_global) {
  return ix_global + Nx[0] * iy_global + Nx[0] * Nx[1] * iz_global;
}