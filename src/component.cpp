// component.cpp
//
// Copyright (c) 2018 Ben Lindsay <benjlindsay@gmail.com>

#include "component.hpp"

Component::Component(Sim *sim) : sim(sim) {
  int M = sim->M;
  rho_center = ArrayXd::Zero(M);
  rho_center_hat = ArrayXd::Zero(M);
  rho = ArrayXd::Zero(M);
  rho_hat = ArrayXd::Zero(M);
}
