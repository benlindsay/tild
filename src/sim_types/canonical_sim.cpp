// canonical_sim.cpp
//
// Copyright (c) 2017 Ben Lindsay <benjlindsay@gmail.com>

#include "canonical_sim.hpp"

Canonical_Sim::Canonical_Sim(YAML::Node input) : Sim(input) {
  utils::print_one_line("Initializing Canonical_Sim");
  init_default_summary_options_list();
  init_output_list(input);
}

void Canonical_Sim::init_default_summary_options_list() {
  default_summary_options_list.push_back("iter");
  default_summary_options_list.push_back("bond_energy");
  default_summary_options_list.push_back("nonbond_energy");
}

void Canonical_Sim::init_output_list(YAML::Node input) {
  fs::path output_dir(input["output_dir"].as<std::string>());
  if (!input["outputs"]) {
    // If no outputs are specified in the input file, just do default outputs
    Output *output =
        new Summary_Output(this, output_dir, default_summary_options_list);
    output_list.push_back(output);
  } else {
    // Otherwise, add an Output object for each output type specified in the
    // YAML file
    YAML::Node outputs_node = input["outputs"];
    for (YAML::const_iterator it = outputs_node.begin();
         it != outputs_node.end(); ++it) {
      std::string output_type = it->first.as<std::string>();
      YAML::Node output_type_params = it->second;  // outputs_node[output_type];
      Output *output = Output_Factory::New_Output(this, input, output_type,
                                                  output_type_params);
      output_list.push_back(output);
    }
  }
}

// Private Functions
