// sim.cpp
//
// Copyright (c) 2017 Ben Lindsay <benjlindsay@gmail.com>

#include "sim.hpp"

void Sim::init(YAML::Node input) {
  utils::print_one_line("Initializing Sim");
  // Initialize dim
  if (!input["dim"]) {
    utils::die("dim keyword not found in input file");
  } else {
    dim = input["dim"].as<int>();
  }

  iter = 0;
  do_fractional_molecules = false;

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

  // Initialize kappa
  if (!input["kappa"]) {
    kappa = -1.0;
  } else {
    kappa = input["kappa"].as<double>();
  }

  // Initialize random variables and objects
  if (!input["random_seed"]) {
    auto now = std::chrono::system_clock::now();
    auto now_ms = std::chrono::time_point_cast<std::chrono::milliseconds>(now);
    random_seed = now_ms.time_since_epoch().count();
  } else {
    random_seed = input["random_seed"].as<unsigned long>();
  }
  std::cout << "Random Seed: " << random_seed << std::endl;
  random_generator.seed(random_seed);
  uniform_dist = std::uniform_real_distribution<double>(0.0, 1.0);
  gaussian_dist = std::normal_distribution<double>(0.0, 1.0);

  // Initialize box/grid variables
  init_box_vars(input);

  // Inizialize conv_function_list with empty arrays since the component
  // constructors will fill this list
  conv_function_list =
      std::vector<ArrayXd>(Component::max_n_species, ArrayXd::Zero(0));

  // Initialize components
  // This will call the derived class's init_component_list() function (i.e.
  // Semi_Grand_Sim::init_component_list() ) instead of the base class
  // ( Sim::init_component_list() ) if one has been defined
  init_component_list(input);

  recalculate_rho_0();

  // Loop over all species in all components to find the number of species in
  // system
  n_species = 0;
  for (size_t i_comp = 0; i_comp < component_list.size(); i_comp++) {
    for (size_t i = 0; i < component_list[i_comp]->species_list.size(); i++) {
      int species = component_list[i_comp]->species_list[i];
      n_species = std::max(n_species, species + 1);
    }
  }

  // Initialize lists that depend on n_species
  conv_function_list.resize(n_species);
  species_list = std::vector<bool>(n_species, 0);
  diffusion_coeff_list = std::vector<double>(n_species, 1.0);
  species_density_list = std::vector<ArrayXd>(n_species, ArrayXd::Zero(ML));
  grad_field_list = std::vector<ArrayXXd>(n_species, ArrayXXd::Zero(ML, dim));

  // Fill species_list
  for (size_t i_comp = 0; i_comp < component_list.size(); i_comp++) {
    for (size_t i = 0; i < component_list[i_comp]->species_list.size(); i++) {
      int species = component_list[i_comp]->species_list[i];
      species_list[species] = true;
    }
  }

  // Check if there are any unused species between 0 and n_species-1
  for (int species = 0; species < n_species; species++) {
    if (species_list[species] == false) {
      std::stringstream ss;
      ss << "Warning: Species up through '"
         << Component::species_int_to_char(n_species - 1)
         << "' are present, but you're missing '"
         << Component::species_int_to_char(species) << "'." << std::endl;
      utils::print_one_line(ss);
    }
  }

  if (input["restart_file"]) {
    fs::path restart_file_path(input["restart_file"].as<std::string>());
    load_restart_file(restart_file_path);
  }

  // Initialize chi
  chi = ArrayXXd::Zero(n_species, n_species);
  if (input["chi"]) {
    for (YAML::const_iterator it = input["chi"].begin();
         it != input["chi"].end(); it++) {
      std::string species_pair = it->first.as<std::string>();
      char species_1_char = species_pair[0];
      char species_2_char = species_pair[1];
      int species_1 = Component::species_char_to_int(species_1_char);
      int species_2 = Component::species_char_to_int(species_2_char);
      // Ensure species_1 < species_2
      if (species_1 > species_2) {
        int tmp = species_1;
        species_1 = species_2;
        species_2 = tmp;
      }
      // Assign chi value to the right spot in the chi array
      chi(species_1, species_2) = it->second.as<double>();
    }
  }
  std::cout << "chi grid:" << std::endl;
  std::cout << chi << std::endl;

  // Read diffusion coefficients from input file
  if (input["diffusion_coeffs"]) {
    for (YAML::const_iterator it = input["diffusion_coeffs"].begin();
         it != input["diffusion_coeffs"].end(); ++it) {
      char species_char = it->first.as<char>();
      double diffusion_coeff = it->second.as<double>();
      int species = Component::species_char_to_int(species_char);
      diffusion_coeff_list[species] = diffusion_coeff;
    }
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
  free(fftw_in_array);
  free(fftw_out_array);
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
        if (Nx[i] % 2 == 0) {
          std::stringstream ss;
          ss << "Nx[" << i << "] was " << Nx[i] << " which is even. "
             << "Adding 1 to make it odd." << std::endl;
          utils::print_one_line(ss);
          Nx[i] += 1;
        }
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
  bool fractional_component = false;
  for (std::size_t i = 0; i < components.size(); i++) {
    Component *comp = Component_Factory::New_Component(
        this, input, components[i], fractional_component);
    component_list.push_back(comp);
  }
}

void Sim::load_restart_file(fs::path restart_file_path) {
  std::string extension = fs::extension(restart_file_path);
  if (extension != ".lammpstrj") {
    utils::die(std::string() + extension +
               " not a recognized restart file type. Only .lammpstrj files "
               "accepted for now.");
  }
  if (!fs::exists(restart_file_path)) {
    utils::die(std::string("Can't find ") + restart_file_path.string() + "!");
  }
  DumpIO dump_in = DumpIO();
  std::ifstream restart_file(restart_file_path.string());
  dump_in.OpenI(restart_file);
  size_t n_frames = dump_in.NumFrames(restart_file);
  for (size_t i = 0; i < n_frames; i++) {
    dump_in.LoadNextFrame(restart_file);
  }
  std::vector<std::vector<double> > positions;
  dump_in.GetPos(positions);
  size_t n_sites = 0;
  for (size_t i = 0; i < component_list.size(); i++) {
    n_sites += component_list[i]->n_sites;
  }
  assert(n_sites == positions.size());
  std::vector<int> mols;
  dump_in.GetDataCol("mol", mols);
  // If mol indexing in restart file starts with something other than 0, shift
  // the mol values to start with 0
  int mol_index_shift = mols[0];
  if (mol_index_shift != 0) {
    for (size_t i = 0; i < mols.size(); i++) {
      mols[i] -= mol_index_shift;
    }
  }
  std::vector<int> types;
  dump_in.GetDataCol("type", types);
  size_t prev_comp_sites = 0;
  size_t prev_comp_molecules = 0;
  for (size_t i_comp = 0; i_comp < component_list.size(); i_comp++) {
    Component *comp = component_list[i_comp];
    for (int i_site = 0; i_site < component_list[i_comp]->n_sites; i_site++) {
      int next_site = prev_comp_sites + i_site;
      int i_mol = comp->get_mol_id(i_site) + prev_comp_molecules;
      // Make sure site id and molecule id match what's expected based on yaml
      // input file
      assert(comp->site_types[i_site] == types[next_site]);
      assert(i_mol == mols[next_site]);
      for (int d = 0; d < dim; d++) {
        comp->site_coords(i_site, d) = positions[next_site][d];
      }
    }
    prev_comp_sites += comp->n_sites;
    prev_comp_molecules += comp->n_molecules;
  }
}

void Sim::init_potentials() {
  for (int species = 0; species < n_species; species++) {
    grad_field_list[species] = ArrayXXd::Zero(ML, dim);
  }

  pair_potential_arrays = std::vector<std::vector<ArrayXd> >(
      n_species, std::vector<ArrayXd>(n_species, ArrayXd()));

  pair_potential_gradient_arrays = std::vector<std::vector<ArrayXXd> >(
      n_species, std::vector<ArrayXXd>(n_species, ArrayXXd()));

  pair_potential_gradient_hat_arrays = std::vector<std::vector<ArrayXXcd> >(
      n_species, std::vector<ArrayXXcd>(n_species, ArrayXXcd()));

  for (int species_i = 0; species_i < n_species; species_i++) {
    ArrayXd conv_func_i = conv_function_list[species_i];
    for (int species_j = species_i; species_j < n_species; species_j++) {
      ArrayXd conv_func_j = conv_function_list[species_j];
      ArrayXd pair_potential = ArrayXd::Zero(ML);
      ArrayXXd pair_potential_gradient = ArrayXXd::Zero(ML, dim);
      ArrayXXcd pair_potential_gradient_hat = ArrayXXcd::Zero(ML, dim);
      convolve(conv_func_i, conv_func_j, pair_potential);
      // Add a factor of V which is necessary for future FFT operations
      pair_potential *= V;
      calculate_gradients(pair_potential, pair_potential_gradient_hat,
                          pair_potential_gradient);

      // Store calculated arrays inside 2d vectors for later access
      pair_potential_arrays[species_i][species_j] = pair_potential;
      pair_potential_gradient_hat_arrays[species_i][species_j] =
          pair_potential_gradient_hat;
      pair_potential_gradient_arrays[species_i][species_j] =
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

void Sim::recalculate_rho_0() {
  double total_component_mass = 0.0;
  for (size_t i_comp = 0; i_comp < component_list.size(); i_comp++) {
    Component *comp = component_list[i_comp];
    double component_mass =
        (comp->n_molecules - 1 + comp->last_molecule_fractional_presence) *
        comp->molecule_mass;
    total_component_mass += component_mass;
  }
  rho_0 = total_component_mass / V;
}

void Sim::recalculate_component_vol_fracs() {
  recalculate_rho_0();
  for (size_t i_comp = 0; i_comp < component_list.size(); i_comp++) {
    Component *comp = component_list[i_comp];
    double component_mass =
        (comp->n_molecules - 1 + comp->last_molecule_fractional_presence) *
        comp->molecule_mass;
    comp->vol_frac = component_mass / (rho_0 * V);
  }
}

void Sim::calculate_grid_densities() {
  // Zero the species density grids stored in species_density_map
  for (int species = 0; species < n_species; species++) {
    // Add each components' rho centers to the total density for each species
    species_density_list[species] = ArrayXd::Zero(ML);
  }
  for (size_t i_comp = 0; i_comp < component_list.size(); i_comp++) {
    Component *comp = component_list[i_comp];
    comp->calculate_grid_densities();
    for (size_t species = 0; species < comp->rho_center_list.size();
         species++) {
      // Add each components' rho centers to the total density for each species
      species_density_list[species] += comp->rho_center_list[species];
    }
  }
}

void Sim::calculate_forces() {
  bond_energy = 0.0;
  nonbond_energy = 0.0;

  // Bonded forces
  for (size_t i_comp = 0; i_comp < component_list.size(); i_comp++) {
    Component *comp = component_list[i_comp];
    // Zero forces array
    comp->site_forces = ArrayXXd::Zero(comp->n_sites, dim);
    bond_energy += comp->calculate_bond_forces_and_energy();
  }

  for (int species = 0; species < n_species; species++) {
    grad_field_list[species] = ArrayXXd::Zero(ML, dim);
  }

  // Fill fields in grad_field_map with species-species nonbonded interactions
  ArrayXcd rho_hat(ML);
  ArrayXXcd field_grad_prod_hat(ML, dim);
  ArrayXXd field_grad_prod(ML, dim);
  std::complex<double> *field_grad_prod_hat_ptr = field_grad_prod_hat.data();
  double *field_grad_prod_ptr = field_grad_prod.data();
  for (int species_i = 0; species_i < n_species; species_i++) {
    fftw_fwd(species_density_list[species_i], rho_hat);
    for (int species_j = 0; species_j < n_species; species_j++) {
      int species_1 = std::min(species_i, species_j);
      int species_2 = std::max(species_i, species_j);
      double factor = 0.0;
      factor += chi(species_1, species_2) / rho_0;
      factor += kappa / rho_0;
      if (factor == 0.0) {
        continue;
      }
      field_grad_prod_hat =
          pair_potential_gradient_hat_arrays[species_1][species_2].colwise() *
          rho_hat;
      for (int d = 0; d < dim; d++) {
        fftw_back(field_grad_prod_hat_ptr + d * ML,
                  field_grad_prod_ptr + d * ML);
      }
      // Species i acting on Species j
      grad_field_list[species_j] += field_grad_prod * factor;
    }
  }

  // Accumulate the nonbonded forces
  for (size_t i_comp = 0; i_comp < component_list.size(); i_comp++) {
    Component *comp = component_list[i_comp];
    for (int i_site = 0; i_site < comp->n_sites; i_site++) {
      int species = comp->site_types[i_site];
      for (int d = 0; d < dim; d++) {
        double force = 0.0;
        for (int i_grid = 0; i_grid < n_subgrid_points; i_grid++) {
          int grid_ind = comp->site_grid_indices(i_site, i_grid);
          double grid_weight = comp->site_grid_weights(i_site, i_grid);
          force -= grad_field_list[species](grid_ind, d) * grid_weight *
                   grid_point_volume;
        }
        if (i_comp < 2 &&
            comp->site_molecule_ids[i_site] == comp->n_molecules - 1) {
          force *= comp->last_molecule_fractional_presence;
        }
        comp->site_forces(i_site, d) += force;
        if (comp->site_forces(i_site, d) > 10000.0) {
          std::cout << "force = " << comp->site_forces(i_site, d) << " > 100"
                    << std::endl;
        }
      }
    }
    // nonbond_energy += calculate_nonbond_energy();
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

void Sim::get_spline_weights(ArrayXd &dx_from_nearest_grid_point,
                             double *axes_grid_weights) {
  ArrayXd dx_norm = dx_from_nearest_grid_point / dx;
  double *dx_norm_ptr = dx_norm.data();

  double scale = double(M) / V;
  if (mesh_order == 0) {
    // TODO: Is this true? Why the scale when grid weights sum to 1 along each
    // dimension in all other cases?
    axes_grid_weights[0] = 1.0 * scale;
  } else if (mesh_order == 1) {
    for (int d = 0; d < dim; d++) {
      int shift = d * (mesh_order + 1);
      double *w = axes_grid_weights + shift;
      double d1 = dx_norm_ptr[d];
      w[0] = (1 - 2 * d1) / 2.0;
      w[1] = (1 + 2 * d1) / 2.0;
    }
  } else if (mesh_order == 2) {
    for (int d = 0; d < dim; d++) {
      int shift = d * (mesh_order + 1);
      double *w = axes_grid_weights + shift;
      double d1 = dx_norm_ptr[d];
      double d2 = d1 * d1;

      w[0] = (1 - 4 * d1 + 4 * d2) / 8.0;
      w[1] = (3 - 4 * d2) / 4.0;
      w[2] = (1 + 4 * d1 + 4 * d2) / 8.0;
    }
  } else if (mesh_order == 3) {
    for (int d = 0; d < dim; d++) {
      int shift = d * (mesh_order + 1);
      double *w = axes_grid_weights + shift;
      double d1 = dx_norm_ptr[d];
      double d2 = d1 * d1;
      double d3 = d2 * d1;

      w[0] = (1 - 6 * d1 + 12 * d2 - 8 * d3) / 48.0;
      w[1] = (23 - 30 * d1 - 12 * d2 + 24 * d3) / 48.0;
      w[2] = (23 + 30 * d1 - 12 * d2 - 24 * d3) / 48.0;
      w[3] = (1 + 6 * d1 + 12 * d2 + 8 * d3) / 48.0;
    }
  } else if (mesh_order == 4) {
    for (int d = 0; d < dim; d++) {
      int shift = d * (mesh_order + 1);
      double *w = axes_grid_weights + shift;
      double d1 = dx_norm_ptr[d];
      double d2 = d1 * d1;
      double d3 = d2 * d1;
      double d4 = d3 * d1;

      w[0] = (1 - 8 * d1 + 24 * d2 - 32 * d3 + 16 * d4) / 384.0;
      w[1] = (19 - 44 * d1 + 24 * d2 + 16 * d3 - 16 * d4) / 96.0;
      w[2] = (115 - 120 * d2 + 48 * d4) / 192.0;
      w[3] = (19 + 44 * d1 + 24 * d2 - 16 * d3 - 16 * d4) / 96.0;
      w[4] = (1 + 8 * d1 + 24 * d2 + 32 * d3 + 16 * d4) / 384.0;
    }
  } else if (mesh_order == 5) {
    for (int d = 0; d < dim; d++) {
      int shift = d * (mesh_order + 1);
      double *w = axes_grid_weights + shift;
      double d1 = dx_norm_ptr[d];
      double d2 = d1 * d1;
      double d3 = d2 * d1;
      double d4 = d3 * d1;
      double d5 = d4 * d1;

      w[0] = (1 - 10 * d1 + 40 * d2 - 80 * d3 + 80 * d4 - 32 * d5) / 3840.0;
      w[1] =
          (237 - 750 * d1 + 840 * d2 - 240 * d3 - 240 * d4 + 160 * d5) / 3840.0;
      w[2] =
          (841 - 770 * d1 - 440 * d2 + 560 * d3 + 80 * d4 - 160 * d5) / 1920.0;
      w[3] =
          (841 + 770 * d1 - 440 * d2 - 560 * d3 + 80 * d4 + 160 * d5) / 1920.0;
      w[4] = (237 + 750 * d1 + 840 * d2 + 2400 * d3 - 240 * d4 - 160 * d5) /
             3840.0;
      w[5] = (1 + 10 * d1 + 40 * d2 + 80 * d3 + 80 * d4 + 32 * d5) / 3840.0;
    }
  } else {
    utils::die("get_spline_weights not set up for this interpolation order!\n");
  }
  for (int i_w = 0; i_w < dim * (mesh_order + 1); i_w++) {
    if (axes_grid_weights[i_w] < 0) {
      utils::die("negative weights");
    }
  }
}
