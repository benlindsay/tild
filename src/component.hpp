// component.hpp
//
// Copyright (c) 2018 Ben Lindsay <benjlindsay@gmail.com>

#ifndef COMPONENT_HPP
#define COMPONENT_HPP

#include <cassert>
#include <fstream>
#include <map>
#include "Eigen/Dense"
#include "globals.hpp"
#include "utils.hpp"

class Sim;

using Eigen::Array;
using Eigen::Dynamic;
using Eigen::RowMajor;
using Eigen::ArrayXi;   // Dynamically sized integer Array
using Eigen::ArrayXd;   // Dynamically sized double Array
using Eigen::ArrayXXd;  // Dynamically sized 2D double Array
typedef Array<double, Dynamic, Dynamic, RowMajor>
    ArrayXXdR;  // Like ArrayXXd but with data contiguous along the rows
typedef Array<int, Dynamic, Dynamic, RowMajor>
    ArrayXXiR;  // Like ArrayXXi but with data contiguous along the rows

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
  static const int max_n_species =
      7;  // Should match number of letters in Species_Type enum above
  void init_site_grid_vars();
  virtual void calculate_axes_grid_weights(int i_site,
                                           ArrayXi &subgrid_center_indices,
                                           ArrayXXdR &axes_grid_weights);
  virtual void add_site_to_grid(int i_site, ArrayXi &subgrid_center_indices,
                                ArrayXXdR &axes_grid_weights);
  virtual void calculate_grid_densities();
  virtual double calculate_bond_forces_and_energy() = 0;
  virtual int get_mol_id(int site_id) = 0;
  virtual void move_particles();

  Sim *sim;
  std::string name;  // Name of component
  double vol_frac;   // Total volume fraction
  int n_molecules;   // Total number of molecules of this type in system
  int n_sites;       // Total number of sites (anything with its own set
                     // of coordinates) of this type in system

  // Store the center density distributions for each species in the component in
  // a map data structure. For example, the B species center density array can
  // be accessed in any of the following ways:
  //   ArrayXd rho_b = rho_center_map[1];
  //   ArrayXd rho_b = rho_center_map[Component::B];
  //   Species_Type b = Component::B; ArrayXd rho_b = rho_center_map[b]
  std::map<Species_Type, ArrayXd> rho_center_map;

  ArrayXi site_types;    // Array of species type for each site
  ArrayXi molecule_ids;  // Each component has ids from 0 to n_molecules so that
                         // they don't have change if other components add or
                         // subtract molecules
  ArrayXXd site_coords;  // Coordinates of sites
  ArrayXXd site_forces;  // Coordinates of sites
  ArrayXXiR site_grid_indices;  // Element (i, j) contains global (i.e. [0, M) )
                                // index of jth grid point near particle i
  ArrayXXdR
      site_grid_weights;  // Element (i, j) contains weight associated with
                          // grid point whos element you see in
                          // site_grid_indices(i, j)
};

#include "component_types/homopolymer.hpp"
#include "sim.hpp"

#endif  // COMPONENT_HPP