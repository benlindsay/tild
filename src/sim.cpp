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

  // Assign V and M
  V = Lx.prod();
  M = Nx.prod();
  ML = M;

  local_grid_coords = ArrayXXd::Zero(ML, dim);
  std::vector<ArrayXd> axes;
  for (int i = 0; i < dim; i++) {
    axes.push_back(ArrayXd::LinSpaced(Nx[i], 0.0, Lx[i] - dx[i]));
  }
  int n_reps = 1;
  for (int d = 0; d < dim; d++) {
    for (int i = 0; i < ML; i++) {
      local_grid_coords(i, d) = axes[d][(i / n_reps) % Nx[d]];
    }
    n_reps *= Nx[d];
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

void Sim::write_iter_0_outputs() {
  for (size_t i = 0; i < output_list.size(); i++) {
    if (output_list[i]->is_time_to_write()) {
      output_list[i]->write_iter_0();
    }
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
  write_iter_0_outputs();
  for (iter = 1; iter <= max_iter; iter++) {
    calculate_grid_densities();
    calculate_forces();
    move_particles();
    write_outputs();
  }
}