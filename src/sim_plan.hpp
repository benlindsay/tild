// sim_plan.hpp
//
// Copyright (c) 2017 Ben Lindsay <benjlindsay@gmail.com>

#ifndef SIM_PLAN_HPP
#define SIM_PLAN_HPP

#include "globals.hpp"
#include "sim.hpp"
#include "sim_factory.hpp"
#include "utils.hpp"
#include "yaml-cpp/yaml.h"

class Sim_Plan {
 public:
  Sim_Plan(YAML::Node input){};
  virtual ~Sim_Plan(){};
  virtual void run(void) = 0;
};

#endif  // SIM_PLAN_HPP
