// canonical_sim.cpp
//
// Copyright (c) 2017 Ben Lindsay <benjlindsay@gmail.com>

#include "canonical_sim.hpp"

Canonical_Sim::Canonical_Sim(YAML::Node input) : Sim(input) {
  utils::print_one_line("Initializing Canonical_Sim");
  init_default_summary_var_list();
  init_component_list(input);
  init_output_list(input);
}

Canonical_Sim::~Canonical_Sim() {
  for (size_t i = 0; i < output_list.size(); i++) {
    delete output_list[i];
  }
}

std::string Canonical_Sim::get_var_as_string(std::string var_name,
                                             int str_len) {
  std::string var_name_lower = utils::to_lower(var_name);
  std::ostringstream os;
  // Assuming we want 1 space between fields, 2 characters go to "#.", and 4
  // characters go to exponent, i.e. "e+00", string length - 7 gives the number
  // of digits we can put after the decimal
  int scientific_precision = str_len - 7;
  // Set number of digits after decimal point in scientific notation
  os.precision(scientific_precision);
  // Set minimum width of output string (left-padded with spaces)
  os << std::setw(str_len - 1);
  // Add scientific notation flag to force specified number of digits after
  // decimal point
  os << std::scientific;
  if (var_name_lower == "iter") {
    os << iter;
  } else if (var_name_lower == "u_bond") {
    os << 42.0;
  } else {
    utils::die("Can't find match for " + var_name);
  }
  return " " + os.str();
}

void Canonical_Sim::init_default_summary_var_list() {
  default_summary_var_list.push_back("iter");
  default_summary_var_list.push_back("u_bond");
}

void Canonical_Sim::init_output_list(YAML::Node input) {
  fs::path output_dir(input["output_dir"].as<std::string>());
  if (!input["outputs"]) {
    // If no outputs are specified in the input file, just do default outputs
    Output *output =
        new Summary_Output(this, output_dir, default_summary_var_list);
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

void Canonical_Sim::calculate_grid_densities() {
  utils::print_one_line("Calculating grid densities");
}

void Canonical_Sim::calculate_forces() {
  utils::print_one_line("Calculating forces");
}

void Canonical_Sim::move_particles() {
  utils::print_one_line("Moving particles");
}

// Private Functions
