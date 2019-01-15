// homopolymer.hpp
//
// Copyright (c) 2018 Ben Lindsay <benjlindsay@gmail.com>

#ifndef HOMOPOLYMER_HPP
#define HOMOPOLYMER_HPP

#include "../component.hpp"

class Homopolymer : public Component {
 public:
  Homopolymer(Sim *sim, double vol_frac, int n_segments_per_molecule,
              int species);
  virtual ~Homopolymer(){};
  virtual double calculate_bond_forces_and_energy();
  virtual int get_mol_id(int site_id);
  int species;
  int n_segments_per_molecule;
};

#endif  // HOMOPOLYMER_HPP
