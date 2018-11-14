// sim_factory.cpp
//
// Copyright (c) 2017 Ben Lindsay <benjlindsay@gmail.com>

#include "sim_factory.hpp"

Sim* Sim_Factory::New_Sim(YAML::Node input) {
  if (!input["sim_type"]) {
    utils::print_one_line(
        "sim_type not included in input."
        "Assuming canonical field theory simulation.");
    return new Canonical_Sim(input);
  } else if (input["sim_type"].as<std::string>() == "canonical") {
    return new Canonical_Sim(input);
  } else {
    utils::die("sim_type " + input["sim_type"].as<std::string>() +
               " not recognized");
    return NULL;
  }
}
