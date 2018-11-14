// brent_plan.hpp
//
// Copyright (c) 2017 Ben Lindsay <benjlindsay@gmail.com>

#ifndef BRENT_PLAN_HPP
#define BRENT_PLAN_HPP

#include <iostream>
#include "../globals.hpp"
#include "../sim.hpp"
#include "../sim_plan.hpp"
#include "yaml-cpp/yaml.h"

class Brent_Plan : public Sim_Plan {
 public:
  Brent_Plan(YAML::Node input);
  virtual ~Brent_Plan() { delete sim; };
  virtual void run(void);
  Sim *sim;
};

#endif  // BRENT_PLAN_HPP
