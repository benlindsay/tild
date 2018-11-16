// homopolymer.hpp
//
// Copyright (c) 2018 Ben Lindsay <benjlindsay@gmail.com>

#ifndef HOMOPOLYMER_HPP
#define HOMOPOLYMER_HPP

#include "../component.hpp"

class Homopolymer : public Component {
 public:
  Homopolymer(Sim *sim, double vol_frac, int n_segments_per_molecule,
              Component::Species_Type species);
  virtual ~Homopolymer(){};
  Component::Species_Type species;
  int n_segments_per_molecule;
};

#endif  // HOMOPOLYMER_HPP