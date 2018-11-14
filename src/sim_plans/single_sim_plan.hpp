// single_sim_plan.hpp
//
// Copyright (c) 2017 Ben Lindsay <benjlindsay@gmail.com>

#ifndef SINGLE_SIM_PLAN_HPP
#define SINGLE_SIM_PLAN_HPP

class Sim_Plan;

#include <iostream>
#include "../sim.hpp"
#include "../sim_plan.hpp"
#include "yaml-cpp/yaml.h"

class Single_Sim_Plan : public Sim_Plan {
 public:
  Single_Sim_Plan(YAML::Node input);
  virtual ~Single_Sim_Plan() { delete sim; };
  virtual void run(void);
  Sim *sim;
};

#endif  // SINGLE_SIM_PLAN_HPP
