// component.cpp
//
// Copyright (c) 2018 Ben Lindsay <benjlindsay@gmail.com>

#include "component.hpp"

Component::Component(Sim *sim, double vol_frac) : sim(sim), vol_frac(vol_frac) {
  int M = sim->M;
  rho_center = ArrayXd::Zero(M);
  rho_center_hat = ArrayXd::Zero(M);
  rho = ArrayXd::Zero(M);
  rho_hat = ArrayXd::Zero(M);
}
