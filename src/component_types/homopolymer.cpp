// ft_homopolymer.cpp
//
// Copyright (c) 2018 Ben Lindsay <benjlindsay@gmail.com>

#include "homopolymer.hpp"

Homopolymer::Homopolymer(Sim *sim, int n_segments_per_molecule,
                         Component::Species_Type species)
    : Component(sim),
      species(species),
      n_segments_per_molecule(n_segments_per_molecule) {
  name = std::string("homopolyer_") + Component::get_species_name(species);
  utils::print_one_line("Initializing Component " + name);
}