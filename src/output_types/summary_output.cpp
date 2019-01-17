// summary_output.cpp
//
// Copyright (c) 2018 Ben Lindsay <benjlindsay@gmail.com>

#include "summary_output.hpp"

Summary_Output::Summary_Output(Sim *sim, std::vector<std::string> options_list,
                               int print_freq, fs::path file_path,
                               int column_width, bool write_header)
    : Output(sim) {
  init(options_list, print_freq, file_path, column_width, write_header);
}

Summary_Output::Summary_Output(Sim *sim, fs::path output_dir,
                               std::vector<std::string> options_list)
    : Output(sim) {
  init(output_dir, options_list);
}

void Summary_Output::init(std::vector<std::string> _options_list,
                          int _print_freq, fs::path _file_path,
                          int _column_width, bool _write_header) {
  options_list = _options_list;
  print_freq = _print_freq;
  file_path = _file_path;
  column_width = _column_width;
  write_header = _write_header;

  // Set width of iteration column to be large enough to handle the last
  // iteration number plus a space on the left side
  iter_column_width = int(std::ceil(std::log10(sim->max_iter))) + 2;
  iter_column_width = std::max(iter_column_width, 5);

  // Create all directories in file_path if they don't exist
  fs::create_directories(file_path.parent_path());
  file.open(file_path);
}

void Summary_Output::init(fs::path output_dir,
                          std::vector<std::string> options_list) {
  int print_freq = Output::default_print_freq;
  fs::path file_path = output_dir / fs::path(Summary_Output::default_file_name);
  int column_width = Output::default_column_width;
  bool write_header = true;

  init(options_list, print_freq, file_path, column_width, write_header);
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
    for (size_t i_opt = 0; i_opt < options_list.size(); i_opt++) {
      std::string opt = options_list[i_opt];
      if (opt == "iter") {
        ss << " " << std::setw(iter_column_width - 1) << opt;
      } else if (opt == "bond_energy") {
        ss << " " << std::setw(column_width - 1) << opt;
      } else if (opt == "nonbond_energy") {
        ss << " " << std::setw(column_width - 1) << opt;
      } else if (opt == "component_vol_fracs") {
        for (size_t i_comp = 0; i_comp < sim->component_list.size(); i_comp++) {
          Component *comp = sim->component_list[i_comp];
          std::string header = std::string("phi_") + comp->abbrev;
          ss << " " << std::setw(column_width - 1) << header;
        }
      } else {
        utils::die("Can't find match for " + opt);
      }
    }
    utils::print_one_line(file, ss);
  }

  // Clear stringstream
  ss.str(std::string());

  for (size_t i_opt = 0; i_opt < options_list.size(); i_opt++) {
    std::string opt = options_list[i_opt];
    if (opt == "iter") {
      ss << " " << std::setw(iter_column_width - 1) << sim->iter;
    } else if (opt == "bond_energy") {
      ss << " ";
      write_double(ss, sim->bond_energy, column_width - 1);
    } else if (opt == "nonbond_energy") {
      ss << " ";
      write_double(ss, sim->nonbond_energy, column_width - 1);
    } else if (opt == "component_vol_fracs") {
      for (size_t i_comp = 0; i_comp < sim->component_list.size(); i_comp++) {
        Component *comp = sim->component_list[i_comp];
        ss << " ";
        write_double(ss, comp->vol_frac, column_width - 1);
      }
    } else {
      utils::die("Can't find match for " + opt);
    }
  }
  utils::print_one_line(file, ss);
}
