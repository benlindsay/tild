// component.hpp
//
// Copyright (c) 2018 Ben Lindsay <benjlindsay@gmail.com>

#ifndef COMPONENT_HPP
#define COMPONENT_HPP

#include "Eigen/Dense"
#include "globals.hpp"

class Sim;

using Eigen::ArrayXi;   // Dynamically sized integer Array
using Eigen::ArrayXd;   // Dynamically sized double Array
using Eigen::ArrayXXd;  // Dynamically sized 2D double Array

class Component {
 public:
  Component(Sim *sim, double vol_frac);
  virtual ~Component(){};
  enum Species_Type { A, B, C, D, E, F, G };
  static char get_species_name(Species_Type s) {
    const char Species_Char[] = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'};
    return Species_Char[s];
  };
  Sim *sim;
  std::string name;  // Name of component
  double vol_frac;   // Total volume fraction
  int n_molecules;   // Total number of molecules of this type in system
  int n_sites;       // Total number of sites (anything with its own set
                     // of coordinates) of this type in system

  ArrayXd rho_center;      // Center density distribution
  ArrayXd rho_center_hat;  // Fourier transform of center density
  ArrayXd rho;             // Smeared/total density distribution
  ArrayXd rho_hat;         // Fourier transform of total density
  ArrayXi site_types;      // Array of species type for each site
  ArrayXi molecule_ids;  // Each component has ids from 0 to n_molecules so that
                         // they don't have change if other components add or
                         // subtract molecules
  ArrayXXd site_coords;  // Coordinates of sites
};

#include "component_types/homopolymer.hpp"
#include "sim.hpp"

#endif  // COMPONENT_HPP