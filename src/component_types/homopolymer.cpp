// homopolymer.cpp
//
// Copyright (c) 2018 Ben Lindsay <benjlindsay@gmail.com>

#include "homopolymer.hpp"

Homopolymer::Homopolymer(Sim *sim, double vol_frac, int n_segments_per_molecule,
                         Component::Species_Type species)
    : Component(sim, vol_frac),
      species(species),
      n_segments_per_molecule(n_segments_per_molecule) {
  name = std::string("homopolyer_") + Component::get_species_name(species);
  utils::print_one_line("Initializing Component " + name);

  n_molecules =
      int(sim->rho_0 * sim->V * vol_frac / n_segments_per_molecule + 0.5);
  n_sites = n_molecules * n_segments_per_molecule;
  site_types = ArrayXi::Constant(n_sites, species);
  molecule_ids = ArrayXi::Zero(n_sites);
  site_coords = ArrayXXd::Zero(n_sites, sim->dim);
  for (int i_mol = 0; i_mol < n_molecules; i_mol++) {
    int i_site_start = i_mol * n_segments_per_molecule;
    molecule_ids.segment(i_site_start, n_segments_per_molecule) = i_mol;
    site_coords.row(i_site_start) = ArrayXXd::Random(1, sim->dim);
    for (int i_seg = 1; i_seg < n_segments_per_molecule; i_seg++) {
      int cur_site = i_site_start + i_seg;
      int prev_site = cur_site - 1;
      site_coords.row(cur_site) =
          site_coords.row(prev_site) + ArrayXXd::Random(1, sim->dim);
    }
  }
}