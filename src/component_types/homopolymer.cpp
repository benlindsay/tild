// homopolymer.cpp
//
// Copyright (c) 2018 Ben Lindsay <benjlindsay@gmail.com>

#include "homopolymer.hpp"

Homopolymer::Homopolymer(Sim *sim, double vol_frac, int n_segments_per_molecule,
                         Component::Species_Type species)
    : Component(sim, vol_frac),
      species(species),
      n_segments_per_molecule(n_segments_per_molecule) {
  name = std::string("homopolyer_") + Component::species_enum_to_char(species);
  utils::print_one_line("Initializing Component " + name);

  // Initialize center density array for specified species to array or zeros
  rho_center_map[species] = ArrayXd::Zero(sim->ML);

  n_molecules =
      int(sim->rho_0 * sim->V * vol_frac / n_segments_per_molecule + 0.5);
  n_sites = n_molecules * n_segments_per_molecule;
  site_types = ArrayXi::Constant(n_sites, species);
  molecule_ids = ArrayXi::Zero(n_sites);
  site_coords = ArrayXXd::Zero(n_sites, sim->dim);
  for (int i_mol = 0; i_mol < n_molecules; i_mol++) {
    int i_site_start = i_mol * n_segments_per_molecule;
    molecule_ids.segment(i_site_start, n_segments_per_molecule) = i_mol;
    site_coords.row(i_site_start) =
        (ArrayXXd::Random(1, sim->dim) + 1.0) * 0.5 * sim->Lx.transpose();
    for (int i_seg = 1; i_seg < n_segments_per_molecule; i_seg++) {
      int cur_site = i_site_start + i_seg;
      int prev_site = cur_site - 1;
      site_coords.row(cur_site) =
          site_coords.row(prev_site) + ArrayXXd::Random(1, sim->dim);
      utils::enforce_pbc(site_coords, sim->Lx, cur_site);
    }
  }
  bool pbc_check = utils::check_pbc(site_coords, sim->Lx);
  std::cout << "PBC Check: " << pbc_check << std::endl;

  double monomer_size_squared = sim->monomer_size * sim->monomer_size;
  double prefactor =
      1.0 / std::pow(2.0 * PI * monomer_size_squared, double(sim->dim) / 2.0);
  ArrayXXd pbc_dist_from_origin = sim->local_grid_coords;
  for (int d = 0; d < sim->dim; d++) {
    for (int i = 0; i < sim->ML; i++) {
      double Lx_d = sim->Lx[d];
      if (pbc_dist_from_origin(i, d) > Lx_d / 2.0) {
        pbc_dist_from_origin(i, d) -= Lx_d;
      }
    }
  }
  ArrayXd dist_from_origin_squared =
      pbc_dist_from_origin.square().rowwise().sum();
  ArrayXd conv_function = prefactor * Eigen::exp(-dist_from_origin_squared /
                                                 (2.0 * monomer_size_squared));
  if (sim->conv_function_map.count(species) > 0) {
    if (conv_function.isApprox(sim->conv_function_map[species])) {
      std::cout << "Species " << Component::species_enum_to_char(species)
                << " conv_function matches one already seen. Good."
                << std::endl;
    } else {
      std::ostringstream s;
      s << "Species " << Component::species_enum_to_char(species)
        << " conv_function doesn't match one already seen. "
        << "Fix that or give it a new species type.";
      utils::die(s.str());
    }
  } else {
    sim->conv_function_map[species] = conv_function;
  }
}

double Homopolymer::calculate_bond_forces_and_energy() {
  double bond_energy = 0.0;
  for (int i_mol = 0; i_mol < n_molecules; i_mol++) {
    for (int i_seg = 0; i_seg < n_segments_per_molecule - 1; i_seg++) {
      int site_1 = i_mol * n_segments_per_molecule + i_seg;
      int site_2 = site_1 + 1;
      ArrayXd dr = sim->pbc_r2_minus_r1(site_coords.row(site_1),
                                        site_coords.row(site_2));
      bond_energy += 1.5 * dr.square().sum();
      site_forces.row(site_1) += 3.0 * dr;
      site_forces.row(site_2) -= 3.0 * dr;
    }
  }
  return bond_energy;
}