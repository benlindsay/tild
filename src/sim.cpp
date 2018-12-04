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

  // Initialize bond_length
  if (!input["bond_length"]) {
    bond_length = 1.0;
  } else {
    bond_length = input["bond_length"].as<double>();
  }

  // Initialize monomer_size
  if (!input["monomer_size"]) {
    monomer_size = 1.0;
  } else {
    monomer_size = input["monomer_size"].as<double>();
  }

  // Initialize timestep
  if (!input["timestep"]) {
    timestep = 1.0;
  } else {
    timestep = input["timestep"].as<double>();
  }

  // Initialize box/grid variables
  init_box_vars(input);

  // Initialize components
  init_component_list(input);

  // Initialize diffusion coefficients
  if (input["diffusion_coeffs"]) {
    for (YAML::const_iterator it = input["diffusion_coeffs"].begin();
         it != input["diffusion_coeffs"].end(); ++it) {
      char species_char = it->first.as<char>();
      double diffusion_coeff = it->second.as<double>();
      Component::Species_Type species =
          Component::species_char_to_enum(species_char);
      diffusion_coeff_map[species] = diffusion_coeff;
    }
  }
  // Set diffusion coefficients to 1.0 if not provided
  for (auto it = conv_function_map.begin(); it != conv_function_map.end();
       it++) {
    Component::Species_Type species = it->first;
    if (diffusion_coeff_map.count(species) == 0) {
      diffusion_coeff_map[species] = 1.0;
    }
  }

  // Set species density grids to zeros
  for (auto it = conv_function_map.begin(); it != conv_function_map.end();
       it++) {
    Component::Species_Type species = it->first;
    species_density_map[species] = ArrayXd::Zero(ML);
  }

  init_potentials();

  // Initialize forces
}

Sim::~Sim() {
  for (size_t i = 0; i < output_list.size(); i++) {
    delete output_list[i];
  }
  for (size_t i = 0; i < component_list.size(); i++) {
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

  init_fftw();

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

  local_grid_k_coords = ArrayXXd::Zero(ML, dim);
  for (int d = 0; d < dim; d++) {
    for (int i = 0; i < ML; i++) {
      int index_i_d = local_grid_indices(i, d);
      if (double(index_i_d) < double(Nx[d]) / 2.0) {
        local_grid_k_coords(i, d) = 2 * PI * double(index_i_d) / Lx[d];
      } else {
        local_grid_k_coords(i, d) = 2 * PI * double(index_i_d - Nx[d]) / Lx[d];
      }
    }
  }

  local_grid_k_magnitude_squared = local_grid_k_coords.square().rowwise().sum();

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

void Sim::init_potentials() {
  for (auto it = conv_function_map.begin(); it != conv_function_map.end();
       it++) {
    Component::Species_Type species = it->first;
    grad_field_map[species] = ArrayXXcd::Zero(ML, dim);
  }

  pair_potential_arrays = std::vector<std::vector<ArrayXd> >(
      Component::max_n_species,
      std::vector<ArrayXd>(Component::max_n_species, ArrayXd()));

  pair_potential_gradient_arrays = std::vector<std::vector<ArrayXXd> >(
      Component::max_n_species,
      std::vector<ArrayXXd>(Component::max_n_species, ArrayXXd()));

  pair_potential_gradient_hat_arrays = std::vector<std::vector<ArrayXXcd> >(
      Component::max_n_species,
      std::vector<ArrayXXcd>(Component::max_n_species, ArrayXXcd()));

  for (auto it_1 = conv_function_map.begin(); it_1 != conv_function_map.end();
       it_1++) {
    int species_int_i = static_cast<int>(it_1->first);
    ArrayXd conv_func_i = it_1->second;
    for (auto it_2 = it_1; it_2 != conv_function_map.end(); it_2++) {
      int species_int_j = static_cast<int>(it_2->first);
      ArrayXd conv_func_j = it_2->second;
      ArrayXd pair_potential = ArrayXd::Zero(ML);
      ArrayXXd pair_potential_gradient = ArrayXXd::Zero(ML, dim);
      ArrayXXcd pair_potential_gradient_hat = ArrayXXcd::Zero(ML, dim);
      convolve(conv_func_i, conv_func_j, pair_potential);
      // Add a factor of V which is necessary for future FFT operations
      pair_potential *= V;
      calculate_gradients(pair_potential, pair_potential_gradient_hat,
                          pair_potential_gradient);

      // Store calculated arrays inside 2d vectors for later access
      pair_potential_arrays[species_int_i][species_int_j] = pair_potential;
      pair_potential_gradient_hat_arrays[species_int_i][species_int_j] =
          pair_potential_gradient_hat;
      pair_potential_gradient_arrays[species_int_i][species_int_j] =
          pair_potential_gradient;
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
  // Write outputs for initial state
  for (iter = 0; iter < max_iter; iter++) {
    calculate_grid_densities();
    calculate_forces();
    write_outputs();
    move_particles();
  }
  write_outputs();
}

void Sim::calculate_grid_densities() {
  for (size_t i_comp = 0; i_comp < component_list.size(); i_comp++) {
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

void Sim::calculate_forces() {
  bond_energy = 0.0;
  nonbond_energy = 0.0;
  for (size_t i_comp = 0; i_comp < component_list.size(); i_comp++) {
    Component *comp = component_list[i_comp];
    // Zero forces array
    comp->site_forces = ArrayXXd::Zero(comp->n_sites, dim);
    bond_energy += comp->calculate_bond_forces_and_energy();
    for (auto it = comp->rho_center_map.begin();
         it != comp->rho_center_map.end(); it++) {
      // it->second points to the actual rho_center density array, for example
      // rho_center_map[A]
      it->second = ArrayXd::Zero(ML);
    }
    comp->calculate_grid_densities();
  }
}

void Sim::move_particles() {
  for (size_t i_comp = 0; i_comp < component_list.size(); i_comp++) {
    Component *comp = component_list[i_comp];
    comp->move_particles();
  }
}

int Sim::get_global_index(int ix_global, int iy_global) {
  return ix_global + Nx[0] * iy_global;
}

int Sim::get_global_index(int ix_global, int iy_global, int iz_global) {
  return ix_global + Nx[0] * iy_global + Nx[0] * Nx[1] * iz_global;
}

ArrayXd Sim::pbc_r2_minus_r1(ArrayXd r1, ArrayXd r2) {
  ArrayXd diff = r2 - r1;
  for (int d = 0; d < dim; d++) {
    if (diff[d] >= Lx[d] / 2.0) {
      diff[d] -= Lx[d];
    } else if (diff[d] < -Lx[d] / 2.0) {
      diff[d] += Lx[d];
    }
  }
  return diff;
}

void Sim::init_fftw() {
  ArrayXi Nx_reverse = Nx.reverse();
  int *Nx_reverse_ptr = Nx_reverse.data();

  ML = M;

  fftw_in_array = (fftw_complex *)fftw_malloc(ML * sizeof(fftw_complex));
  fftw_out_array = (fftw_complex *)fftw_malloc(ML * sizeof(fftw_complex));

  forward_plan = fftw_plan_dft(dim, Nx_reverse_ptr, fftw_in_array,
                               fftw_out_array, FFTW_FORWARD, FFTW_MEASURE);
  backward_plan = fftw_plan_dft(dim, Nx_reverse_ptr, fftw_in_array,
                                fftw_out_array, FFTW_BACKWARD, FFTW_MEASURE);
}

void Sim::fftw_fwd(ArrayXd &in_array, ArrayXcd &out_array) {
  double *in_array_ptr = in_array.data();
  std::complex<double> *out_array_ptr = out_array.data();
  for (int i = 0; i < ML; i++) {
    fftw_in_array[i][0] = in_array_ptr[i];
    fftw_in_array[i][1] = 0.0;
  }
  fftw_execute(forward_plan);
  double norm = 1.0 / double(M);
  for (int i = 0; i < ML; i++) {
    out_array_ptr[i] =
        std::complex<double>(fftw_out_array[i][0], fftw_out_array[i][1]) * norm;
  }
}

void Sim::fftw_back(ArrayXcd &in_array, ArrayXd &out_array) {
  std::complex<double> *in_array_ptr = in_array.data();
  double *out_array_ptr = out_array.data();
  for (int i = 0; i < ML; i++) {
    fftw_in_array[i][0] = in_array_ptr[i].real();
    fftw_in_array[i][1] = in_array_ptr[i].imag();
  }
  fftw_execute(backward_plan);
  for (int i = 0; i < ML; i++) {
    out_array_ptr[i] = fftw_out_array[i][0];
  }
}

void Sim::fftw_back(std::complex<double> *in_array_ptr, double *out_array_ptr) {
  for (int i = 0; i < ML; i++) {
    fftw_in_array[i][0] = in_array_ptr[i].real();
    fftw_in_array[i][1] = in_array_ptr[i].imag();
  }
  fftw_execute(backward_plan);
  for (int i = 0; i < ML; i++) {
    out_array_ptr[i] = fftw_out_array[i][0];
  }
}

void Sim::convolve(ArrayXd &array_1, ArrayXd &array_2, ArrayXd &convolved) {
  ArrayXcd array_1_hat(ML);
  ArrayXcd array_2_hat(ML);
  ArrayXcd convolved_hat(ML);
  fftw_fwd(array_1, array_1_hat);
  fftw_fwd(array_2, array_2_hat);
  convolved_hat = array_1_hat * array_2_hat * V;
  fftw_back(convolved_hat, convolved);
}

void Sim::calculate_gradients(ArrayXd &array, ArrayXXd &grad_arrays) {
  ArrayXcd array_hat(ML);
  fftw_fwd(array, array_hat);
  ArrayXXcd grad_hat_arrays =
      // I * local_grid_k_coords.colwise() *
      local_grid_k_coords.cast<std::complex<double> >().colwise() * array_hat *
      I;
  std::complex<double> *grad_hat_arrays_ptr = grad_hat_arrays.data();
  double *grad_arrays_ptr = grad_arrays.data();
  for (int d = 0; d < dim; d++) {
    // Do the backwards transformation one dimension at a time
    std::complex<double> *in_array_ptr = grad_hat_arrays_ptr + d * ML;
    double *out_array_ptr = grad_arrays_ptr + d * ML;
    fftw_back(in_array_ptr, out_array_ptr);
  }
}

void Sim::calculate_gradients(ArrayXd &array, ArrayXXcd &grad_hat_arrays,
                              ArrayXXd &grad_arrays) {
  ArrayXcd array_hat(ML);
  fftw_fwd(array, array_hat);
  grad_hat_arrays =
      local_grid_k_coords.cast<std::complex<double> >().colwise() * array_hat *
      I;
  std::complex<double> *grad_hat_arrays_ptr = grad_hat_arrays.data();
  double *grad_arrays_ptr = grad_arrays.data();
  for (int d = 0; d < dim; d++) {
    // Do the backwards transformation one dimension at a time
    std::complex<double> *in_array_ptr = grad_hat_arrays_ptr + d * ML;
    double *out_array_ptr = grad_arrays_ptr + d * ML;
    fftw_back(in_array_ptr, out_array_ptr);
  }
}