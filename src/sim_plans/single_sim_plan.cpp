// single_sim_plan.cpp
//
// Copyright (c) 2017 Ben Lindsay <benjlindsay@gmail.com>

#include "single_sim_plan.hpp"

Single_Sim_Plan::Single_Sim_Plan(YAML::Node input) : Sim_Plan(input) {
  utils::print_one_line("Initializing Single_Sim_Plan");
  sim = Sim_Factory::New_Sim(input);
}

void Single_Sim_Plan::run() {
  utils::print_one_line("Running Single_Sim_Plan");
  sim->run();
}
