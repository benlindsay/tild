// component.hpp
//
// Copyright (c) 2018 Ben Lindsay <benjlindsay@gmail.com>

#ifndef COMPONENT_HPP
#define COMPONENT_HPP

#include <unordered_map>
#include "Eigen/Dense"
#include "globals.hpp"
#include "utils.hpp"

class Sim;

using Eigen::ArrayXi;   // Dynamically sized integer Array
using Eigen::ArrayXd;   // Dynamically sized double Array
using Eigen::ArrayXXd;  // Dynamically sized 2D double Array

class Component {
 public:
  Component(Sim *sim, double vol_frac) : sim(sim), vol_frac(vol_frac){};
  virtual ~Component(){};
  enum Species_Type { A, B, C, D, E, F, G };
  static char species_enum_to_char(Species_Type s) {
    const char Species_Char[] = {'a', 'b', 'c', 'd', 'e', 'f', 'g'};
    return Species_Char[s];
  };
  static Species_Type species_char_to_enum(char c) {
    if (int(c) < int('a') || int(c) > int('g')) {
      utils::die(std::string(1, c) + " is not a recognized Species_Type");
    }
    Species_Type species = static_cast<Species_Type>(int(c) - int('a'));
    return species;
  };
  Sim *sim;
  std::string name;  // Name of component
  double vol_frac;   // Total volume fraction
  int n_molecules;   // Total number of molecules of this type in system
  int n_sites;       // Total number of sites (anything with its own set
                     // of coordinates) of this type in system

  // Store the center density distributions for each species in the component in
  // an unordered_map. For example, the B species center density array can be
  // accessed in any of the following ways:
  //   ArrayXd rho_b = rho_center[1];
  //   ArrayXd rho_b = rho_center[Component::B];
  //   Species_Type b = Component::B; ArrayXd rho_b = rho_center[b]
  // The hash<int> thing makes this work with an enum type.
  // See https://stackoverflow.com/a/22361511/2680824
  std::unordered_map<Species_Type, ArrayXd, std::hash<int> > rho_center;
  ArrayXi site_types;    // Array of species type for each site
  ArrayXi molecule_ids;  // Each component has ids from 0 to n_molecules so that
                         // they don't have change if other components add or
                         // subtract molecules
  ArrayXXd site_coords;  // Coordinates of sites
};

#include "component_types/homopolymer.hpp"
#include "sim.hpp"

#endif  // COMPONENT_HPP