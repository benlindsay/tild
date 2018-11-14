// component.hpp
//
// Copyright (c) 2018 Ben Lindsay <benjlindsay@gmail.com>

#ifndef COMPONENT_HPP
#define COMPONENT_HPP

#include "Eigen/Dense"
#include "globals.hpp"

class Sim;

using Eigen::ArrayXcd;  // Dynamically sized complex double Array
using Eigen::ArrayXd;   // Dynamically sized double Array

class Component {
 public:
  Component(Sim *sim) : sim(sim){};
  virtual ~Component(){};
  enum Species_Type { A, B, C, D, E, F, G };
  Sim *sim;
  double vol_frac;         // Total volume fraction
  ArrayXd rho_center;      // Center density distribution
  ArrayXd rho_center_hat;  // Fourier transform of center density
  ArrayXd rho;             // Smeared/total density distribution
  ArrayXd rho_hat;         // Fourier transform of total density
};

#include "sim.hpp"
// #include "component_types/ft_homopolymer.hpp"

#endif  // COMPONENT_HPP