// lammpstrj_output.cpp
//
// Copyright (c) 2018 Ben Lindsay <benjlindsay@gmail.com>

#include "lammpstrj_output.hpp"

Lammpstrj_Output::Lammpstrj_Output(Sim *sim,
                                   std::vector<Component *> component_list,
                                   fs::path output_dir, int print_freq,
                                   std::string name, bool one_frame_per_file)
    : Output(sim) {
  init(component_list, output_dir, print_freq, name, one_frame_per_file);
}

Lammpstrj_Output::~Lammpstrj_Output() { file.close(); }

void Lammpstrj_Output::init(std::vector<Component *> _component_list,
                            fs::path _output_dir, int _print_freq,
                            std::string _name, bool _one_frame_per_file) {
  component_list = _component_list;
  print_freq = _print_freq;
  name = _name;
  one_frame_per_file = _one_frame_per_file;

  file_path = _output_dir / (name + ".lammpstrj");
  file.open(file_path);
}

const std::string Lammpstrj_Output::default_name = "trajectory";

bool Lammpstrj_Output::is_time_to_write(void) {
  if (sim->iter % print_freq == 0 || sim->iter == sim->max_iter) {
    return true;
  } else {
    return false;
  }
}

void Lammpstrj_Output::write_iter_0() {
  // std::stringstream ss;
  // ss << "Hello lammpstrj!";
  // utils::print_one_line(file, ss);
  write();
}

void Lammpstrj_Output::write() {
  std::stringstream ss;
  ss << "iter " << sim->iter;
  utils::print_one_line(file, ss);
}