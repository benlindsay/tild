// lammpstrj_output.cpp
//
// Copyright (c) 2018 Ben Lindsay <benjlindsay@gmail.com>

#include "lammpstrj_output.hpp"

Lammpstrj_Output::Lammpstrj_Output(Sim *sim,
                                   std::vector<Component *> component_list,
                                   fs::path output_dir, int print_freq,
                                   std::string name, bool one_frame_per_file)
    : Output(sim) {
  init(component_list, output_dir, print_freq, name, one_frame_per_file);
}

Lammpstrj_Output::~Lammpstrj_Output() { file.close(); }

void Lammpstrj_Output::init(std::vector<Component *> _component_list,
                            fs::path _output_dir, int _print_freq,
                            std::string _name, bool _one_frame_per_file) {
  component_list = _component_list;
  print_freq = _print_freq;
  name = _name;
  one_frame_per_file = _one_frame_per_file;

  file_path = _output_dir / (name + ".lammpstrj");
  file.open(file_path);
}

const std::string Lammpstrj_Output::default_name = "trajectory";

bool Lammpstrj_Output::is_time_to_write(void) {
  if (sim->iter % print_freq == 0 || sim->iter == sim->max_iter) {
    return true;
  } else {
    return false;
  }
}

void Lammpstrj_Output::write() {
  int n_sites_total = 0;
  for (size_t i = 0; i < component_list.size(); i++) {
    n_sites_total += component_list[i]->n_sites;
  }

  if (RANK == 0) {
    file << "ITEM: TIMESTEP\n";
    file << sim->iter << "\n";
    file << "ITEM: NUMBER OF ATOMS\n";
    file << n_sites_total << "\n";
    file << "ITEM: BOX BOUNDS\n";
    for (int d = 0; d < sim->dim; d++) {
      file << 0.0 << " " << sim->Lx[d] << "\n";
    }
    // lammpstrj files need to be 3D so add a dummy 3rd dimension if simulation
    // is in 2D
    if (sim->dim == 2) {
      file << 0.0 << " " << 1.0 << "\n";
    }
    file << "ITEM: ATOMS id type mol x y z\n";
    int site_id_shift = 0;
    int mol_id_shift = 0;
    for (size_t i_comp = 0; i_comp < component_list.size(); i_comp++) {
      Component *comp = component_list[i_comp];
      for (int i_site = 0; i_site < comp->n_sites; i_site++) {
        file << i_site + site_id_shift << " ";
        file << comp->site_types[i_site] << " ";
        file << comp->molecule_ids[i_site] + mol_id_shift << " ";
        file << comp->site_coords.row(i_site);
        // Dummy 3rd dimension if simulation is 2D
        if (sim->dim == 2) {
          file << " 0.0";
        }
        file << "\n";
      }
      site_id_shift += comp->n_sites;
      mol_id_shift += comp->n_molecules;
    }
  }
}