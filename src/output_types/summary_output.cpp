// summary_output.cpp
//
// Copyright (c) 2018 Ben Lindsay <benjlindsay@gmail.com>

#include "summary_output.hpp"

Summary_Output::Summary_Output(Sim *sim, std::vector<std::string> var_list,
                               int print_freq, fs::path file_path,
                               int column_width, bool write_header)
    : Output(sim) {
  init(sim, var_list, print_freq, file_path, column_width, write_header);
}

Summary_Output::Summary_Output(Sim *sim, fs::path output_dir,
                               std::vector<std::string> var_list)
    : Output(sim) {
  init(sim, output_dir, var_list);
}

void Summary_Output::init(Sim *_sim, std::vector<std::string> _var_list,
                          int _print_freq, fs::path _file_path,
                          int _column_width, bool _write_header) {
  var_list = _var_list;
  print_freq = _print_freq;
  file_path = _file_path;
  column_width = _column_width;
  write_header = _write_header;

  // Create all directories in file_path if they don't exist
  fs::create_directories(file_path.parent_path());
  file.open(file_path);
}

void Summary_Output::init(Sim *sim, fs::path output_dir,
                          std::vector<std::string> var_list) {
  int print_freq = Output::default_print_freq;
  fs::path file_path = output_dir / fs::path(Summary_Output::default_file_name);
  int column_width = Output::default_column_width;
  bool write_header = true;

  init(sim, var_list, print_freq, file_path, column_width, write_header);
}

Summary_Output::~Summary_Output() { file.close(); }

const std::string Summary_Output::default_file_name = "summary.dat";

bool Summary_Output::is_time_to_write() {
  if (sim->iter % print_freq == 0 || sim->iter == sim->max_iter) {
    return true;
  } else {
    return false;
  }
}

void Summary_Output::write() {
  std::stringstream ss;
  if (write_header && sim->iter == 0) {
    for (size_t i = 0; i < var_list.size(); i++) {
      ss << " " << std::setw(column_width - 1) << var_list[i];
    }
    utils::print_one_line(file, ss);
  }
  // Clear stringstream
  ss.str(std::string());
  for (size_t i = 0; i < var_list.size(); i++) {
    ss << sim->get_var_as_string(var_list[i], column_width);
  }
  utils::print_one_line(file, ss);
}
