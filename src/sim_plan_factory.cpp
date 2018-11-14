// sim_plan_factory.cpp
//
// Copyright (c) 2017 Ben Lindsay <benjlindsay@gmail.com>

#include "sim_plan_factory.hpp"

Sim_Plan* Sim_Plan_Factory::New_Sim_Plan(YAML::Node input) {
  if (input["output_dir"]) {
    // If "output_dir" is provided in input file, convert that to an absolute
    // path if it isn't already. If it's a relative path, treat that as relative
    // from the input directory (the directory containing the input file)
    fs::path absolute_output_dir =
        fs::absolute(input["output_dir"].as<std::string>(),
                     input["input_dir"].as<std::string>());
    input["output_dir"] = absolute_output_dir.string();
  } else {
    // If "output_dir" is not provided in the input data, add it to the input
    // data.
    input["output_dir"] = fs::current_path().string();
  }

  if (!input["sim_plan"]) {
    utils::print_one_line(
        "sim_plan not included in input. Assuming single_sim.");
    return new Single_Sim_Plan(input);
  }
  std::string sim_plan_type = input["sim_plan"].as<std::string>();
  if (sim_plan_type == "single_sim") {
    return new Single_Sim_Plan(input);
  } else if (sim_plan_type == "brent") {
    return new Brent_Plan(input);
  } else {
    utils::die("sim_plan " + sim_plan_type + " not recognized");
    return NULL;
  }
}

Sim_Plan* Sim_Plan_Factory::New_Sim_Plan(fs::path input_file_path) {
  // Turn the input_file_path into an absolute path with no ..'s
  input_file_path = fs::canonical(fs::absolute(input_file_path));
  // Get the directory of the input file
  fs::path input_file_dir = input_file_path.parent_path();
  // Read input file and store data in YAML Node
  YAML::Node input = YAML::LoadFile(input_file_path.string());
  // Convert all input data to lowercase
  utils::to_lower(input);
  // Whether or not "input_dir" was provided in the yaml input file, add a
  // variable for input directory to the yaml data
  input["input_dir"] = input_file_dir.string();

  if (!input["output_dir"]) {
    // If "output_dir" is not provided in input file, use the input_dir as the
    // output_dir as well
    input["output_dir"] = input_file_dir.string();
  }

  return New_Sim_Plan(input);
}
