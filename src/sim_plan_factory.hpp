// sim_plan_factory.hpp
//
// Copyright (c) 2017 Ben Lindsay <benjlindsay@gmail.com>

#ifndef SIM_PLAN_FACTORY_HPP
#define SIM_PLAN_FACTORY_HPP

#include "sim_plan.hpp"
#include "sim_plans/brent_plan.hpp"
#include "sim_plans/single_sim_plan.hpp"
#include "utils.hpp"
#include "yaml-cpp/yaml.h"

class Sim_Plan_Factory {
 public:
  static Sim_Plan *New_Sim_Plan(YAML::Node input);
  static Sim_Plan *New_Sim_Plan(fs::path input_file_path);
};

#endif  // SIM_PLAN_FACTORY_HPP
