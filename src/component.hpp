// component.hpp
//
// Copyright (c) 2018 Ben Lindsay <benjlindsay@gmail.com>

#ifndef COMPONENT_HPP
#define COMPONENT_HPP

#include "Eigen/Dense"
#include "globals.hpp"

class Sim;

using Eigen::ArrayXd;   // Dynamically sized double Array
using Eigen::ArrayXXd;  // Dynamically sized 2D double Array

class Component {
 public:
  Component(Sim *sim);
  virtual ~Component(){};
  enum Species_Type { A, B, C, D, E, F, G };
  static char get_species_name(Species_Type s) {
    const char Species_Char[] = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'};
    return Species_Char[s];
  };
  std::string name;  // Name of component
  double vol_frac;   // Total volume fraction
  int n_molecules;   // Total number of molecules of this type in system
  int n_sites;       // Total number of sites (anything with its own set
                     // of coordinates) of this type in system

  Sim *sim;
  ArrayXd rho_center;      // Center density distribution
  ArrayXd rho_center_hat;  // Fourier transform of center density
  ArrayXd rho;             // Smeared/total density distribution
  ArrayXd rho_hat;         // Fourier transform of total density
  ArrayXXd site_coords;    // Coordinates of sites
};

#include "component_types/homopolymer.hpp"
#include "sim.hpp"

#endif  // COMPONENT_HPP