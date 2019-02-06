// sim_factory.cpp
//
// Copyright (c) 2017 Ben Lindsay <benjlindsay@gmail.com>

#include "sim_factory.hpp"

Sim* Sim_Factory::New_Sim(YAML::Node input) {
  Sim* sim;
  if (!input["sim_type"]) {
    utils::print_one_line(
        "sim_type not included in input. "
        "Assuming canonical simulation.");
    sim = new Canonical_Sim();
  } else if (input["sim_type"].as<std::string>() == "canonical") {
    sim = new Canonical_Sim();
  } else if (input["sim_type"].as<std::string>() == "semi_grand") {
    sim = new Semi_Grand_Sim();
  } else {
    utils::die("sim_type " + input["sim_type"].as<std::string>() +
               " not recognized");
    return NULL;
  }
  // Initialize simulation with init() function instead of directly in the
  // constructor so we can access child functions from the parent initialization
  // function. https://stackoverflow.com/a/962148/2680824 explains that you
  // can't do this in a constructor.
  sim->init(input);
  return sim;
}
