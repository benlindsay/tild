// semi_grand_sim.cpp
//
// Copyright (c) 2019 Ben Lindsay <benjlindsay@gmail.com>

#include "semi_grand_sim.hpp"

void Semi_Grand_Sim::init(YAML::Node input) {
  Sim::init(input);
  description = "Semi_Grand_Sim";
  utils::print_one_line("Initializing " + description);
  do_fractional_molecules = true;
  init_default_summary_options_list();
  init_output_list(input);

  if (component_list.size() < 2) {
    utils::die(
        "At least 2 components are required for a Semi-Grand simulation.");
  } else {
    // For now assume the first two components are the ones involved in
    // continuous fractional components and the volume fraction of any other
    // components is fixed
    // swappable_components_list.push_back(0);
    // swappable_components_list.push_back(1);
  }

  // process partial step rate
  if (!input["partial_step_rate"]) {
    utils::die(
        "partial_step_rate necessary in input file for Semi-Grand "
        "simulations.");
  } else {
    partial_step_rate = input["partial_step_rate"].as<double>();
  }
}

void Semi_Grand_Sim::init_default_summary_options_list() {
  default_summary_options_list.push_back("iter");
  default_summary_options_list.push_back("bond_energy");
  default_summary_options_list.push_back("nonbond_energy");
  default_summary_options_list.push_back("vol_fracs");
}

void Semi_Grand_Sim::init_output_list(YAML::Node input) {
  fs::path output_dir(input["output_dir"].as<std::string>());
  if (!input["outputs"]) {
    // If no outputs are specified in the input file, just do default outputs
    Output *output =
        new Summary_Output(this, output_dir, default_summary_options_list);
    output_list.push_back(output);
  } else {
    // Otherwise, add an Output object for each output type specified in the
    // YAML file
    YAML::Node outputs_node = input["outputs"];
    for (YAML::const_iterator it = outputs_node.begin();
         it != outputs_node.end(); ++it) {
      std::string output_type = it->first.as<std::string>();
      YAML::Node output_type_params = it->second;  // outputs_node[output_type];
      Output *output = Output_Factory::New_Output(this, input, output_type,
                                                  output_type_params);
      output_list.push_back(output);
    }
  }
}

void Semi_Grand_Sim::init_component_list(YAML::Node input) {
  utils::print_one_line("Running Semi_Grand_Sim::init_component_list()");
  if (component_list.size() > 0) {
    utils::die("component_list has already been filled!");
  }
  YAML::Node components = input["components"];
  bool fractional_component = false;
  for (size_t i_comp = 0; i_comp < components.size(); i_comp++) {
    // First 2 components in input file are treated as swappable fractional
    // components
    if (i_comp < 2) {
      fractional_component = true;
    } else {
      fractional_component = false;
    }
    Component *comp = Component_Factory::New_Component(
        this, input, components[i_comp], fractional_component);
    component_list.push_back(comp);
  }
  Component *comp_1 = component_list[0];
  Component *comp_2 = component_list[1];
  double total_swappable_mass =
      (comp_1->vol_frac + comp_2->vol_frac) * rho_0 * V;
  comp_1->max_n_molecules = total_swappable_mass / comp_1->molecule_mass;
  comp_2->max_n_molecules = total_swappable_mass / comp_2->molecule_mass;
}

double Semi_Grand_Sim::calculate_dU_chi_kappa_dlambda_1_m_1(Component *comp_1,
                                                            Component *comp_2) {
  ArrayXcd species_density_hat(ML);
  ArrayXcd u_conv_species_density_hat(ML);

  std::vector<std::vector<ArrayXd> > u_conv_species_density_arrays(
      n_species, std::vector<ArrayXd>(n_species, ArrayXd()));

  double dU_chi_kappa_dlambda_1_m_1 = 0.0;

  for (int i_comp = 0; i_comp < 2; i_comp++) {
    Component *comp;
    if (i_comp == 0) {
      comp = comp_1;
    } else if (i_comp == 1) {
      comp = comp_2;
    }
    double component_contribution = 0.0;
    int site_1 = comp->get_mol_first_site_id(comp->n_molecules - 1);
    int n_sites_per_molecule = comp->n_sites / comp->n_molecules;
    for (int species_i = 0; species_i < n_species; species_i++) {
      for (int site_j = site_1; site_j < site_1 + n_sites_per_molecule;
           site_j++) {
        int species_j = comp->site_types[site_j];
        int species_1 = std::min(species_i, species_j);
        int species_2 = std::max(species_i, species_j);
        double chi_plus_kappa = chi(species_1, species_2) + kappa;
        ArrayXd &u_conv_species_density =
            u_conv_species_density_arrays[species_1][species_2];
        if (u_conv_species_density.size() == 0) {
          u_conv_species_density = ArrayXd::Zero(ML);
          convolve(species_density_list[species_i],
                   pair_potential_arrays[species_1][species_2],
                   u_conv_species_density);
        }
        ArrayXi site_j_grid_indices = comp->site_grid_indices.row(site_j);
        ArrayXd site_j_grid_weights = comp->site_grid_weights.row(site_j);
        for (int i_grid = 0; i_grid < site_j_grid_indices.size(); i_grid++) {
          int grid_index = site_j_grid_indices[i_grid];
          component_contribution += chi_plus_kappa *
                                    u_conv_species_density[grid_index] *
                                    site_j_grid_weights[i_grid];
        }
      }  // for site_j
    }    // for species_i
    component_contribution *= grid_point_volume / comp->molecule_mass / rho_0;
    if (i_comp == 1) {
      component_contribution *= -1;
    }
    dU_chi_kappa_dlambda_1_m_1 += component_contribution;
  }  // for i_comp
  return dU_chi_kappa_dlambda_1_m_1;
}

void Semi_Grand_Sim::update_fractional_presence() {
  Component *comp_1 = component_list[0];
  Component *comp_2 = component_list[1];
  // double debroglie_wavelength = 1.0;
  double dH_dlambda_1_m_1 = 0.0;
  dH_dlambda_1_m_1 += calculate_dU_chi_kappa_dlambda_1_m_1(comp_1, comp_2);
  dH_dlambda_1_m_1 -= comp_1->chemical_potential / comp_1->molecule_mass;
  dH_dlambda_1_m_1 += comp_2->chemical_potential / comp_2->molecule_mass;
  // dH_dlambda_1_m_1 += (std::log(comp_1->n_molecules +
  //                               comp_1->last_molecule_fractional_presence)
  //                               +
  //                      dim * std::log(debroglie_wavelength)) /
  //                     comp_1->molecule_mass;
  // dH_dlambda_1_m_1 -= (std::log(comp_2->n_molecules +
  //                               comp_2->last_molecule_fractional_presence)
  //                               +
  //                      dim * std::log(debroglie_wavelength)) /
  //                     comp_2->molecule_mass;
  double noise = gaussian_rand() * std::sqrt(2 * partial_step_rate * timestep);
  double dlambda_1_m_1_dt = -partial_step_rate * dH_dlambda_1_m_1 + noise;
  double delta_lambda_1 = dlambda_1_m_1_dt * timestep / comp_1->molecule_mass;
  double delta_lambda_2 = -dlambda_1_m_1_dt * timestep / comp_2->molecule_mass;
  if (std::abs(delta_lambda_1) > 1.0 || std::abs(delta_lambda_2) > 1.0) {
    std::stringstream ss;
    ss << "delta_lambda_1 = " << delta_lambda_1
       << " and delta_lambda_2 = " << delta_lambda_2
       << ". Their magnitudes should be < 1.0. Try reducing "
          "partial_step_rate or timestep.";
    utils::die(ss);
  }
  comp_1->add_to_fractional_presence(delta_lambda_1);
  comp_2->add_to_fractional_presence(delta_lambda_2);

  // Have different conditions based on what type of components 1 and 2 are

  // Things to do when/after changing volume fractions:
  //   - update species_vol_frac_list
}

void Semi_Grand_Sim::run() {
  utils::print_one_line("Running " + description);
  // Write outputs for initial state
  for (iter = 0; iter < max_iter; iter++) {
    calculate_grid_densities();
    calculate_forces();
    update_fractional_presence();
    write_outputs();
    move_particles();
  }
  write_outputs();
}
