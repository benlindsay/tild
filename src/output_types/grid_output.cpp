// grid_output.cpp
//
// Copyright (c) 2018 Ben Lindsay <benjlindsay@gmail.com>

#include "grid_output.hpp"

Grid_Output::Grid_Output(Sim *sim, fs::path output_dir, int print_freq,
                         int column_width, bool write_header,
                         bool pm3d_compatible)
    : Output(sim) {
  init(output_dir, print_freq, column_width, write_header, pm3d_compatible);
}

Grid_Output::Grid_Output(Sim *sim, fs::path output_dir) : Output(sim) {
  init(output_dir);
}

void Grid_Output::init(fs::path _output_dir, int _print_freq, int _column_width,
                       bool _write_header, bool _pm3d_compatible) {
  output_dir = _output_dir;
  print_freq = _print_freq;
  column_width = _column_width;
  write_header = _write_header;
  pm3d_compatible = _pm3d_compatible;
  if (sim->dim == 3) {
    // No reason to put blank lines in the file for 3D data
    pm3d_compatible = false;
  }
}

void Grid_Output::init(fs::path output_dir) {
  int print_freq = Output::default_print_freq;
  int column_width = Output::default_column_width;
  bool write_header = false;
  bool pm3d_compatible = true;
  init(output_dir, print_freq, column_width, write_header, pm3d_compatible);
}

bool Grid_Output::is_time_to_write() {
  if (sim->iter % print_freq == 0) {
    return true;
  } else {
    return false;
  }
}

void Grid_Output::write_one_file(fs::path file_path, ArrayXd &data) {
  file.open(file_path);
  if (write_header && !pm3d_compatible) {
    file << " " << std::setw(column_width - 1) << "x";
    file << " " << std::setw(column_width - 1) << "y";
    if (sim->dim == 3) {
      file << " " << std::setw(column_width - 1) << "z";
    }
    file << " " << std::setw(column_width - 1) << "data"
         << "\n";
  }
  // Set number of digits after decimal point in scientific notation
  // Assuming we want 1 space between fields, 2 characters go to "#.", and 4
  // characters go to exponent, i.e. "e+00", string length - 7 gives the number
  // of digits we can put after the decimal
  file.precision(column_width - 7);
  file << std::scientific;
  for (int i = 0; i < sim->ML; i++) {
    for (int d = 0; d < sim->dim; d++) {
      file << " " << std::setw(column_width - 1)
           << sim->local_grid_coords(i, d);
    }
    file << " " << data[i] << "\n";
    if (pm3d_compatible && (i + 1) % sim->Nx[0] == 0) {
      file << "\n";
    }
  }
  file.close();
}

void Grid_Output::write_one_file(fs::path file_path, ArrayXcd &data) {
  file.open(file_path);
  if (write_header && !pm3d_compatible) {
    file << " " << std::setw(column_width - 1) << "x";
    file << " " << std::setw(column_width - 1) << "y";
    if (sim->dim == 3) {
      file << " " << std::setw(column_width - 1) << "z";
    }
    file << " " << std::setw(column_width - 1) << "real";
    file << " " << std::setw(column_width - 1) << "imag"
         << "\n";
  }
  // Set number of digits after decimal point in scientific notation
  // Assuming we want 1 space between fields, 2 characters go to "#.", and 4
  // characters go to exponent, i.e. "e+00", string length - 7 gives the number
  // of digits we can put after the decimal
  file.precision(column_width - 7);
  file << std::scientific;
  for (int i = 0; i < sim->ML; i++) {
    for (int d = 0; d < sim->dim; d++) {
      file << " " << std::setw(column_width - 1)
           << sim->local_grid_coords(i, d);
    }
    file << " " << data[i].real();
    file << " " << data[i].imag() << "\n";
    if (pm3d_compatible && (i + 1) % sim->Nx[0] == 0) {
      file << "\n";
    }
  }
  file.close();
}

void Grid_Output::write() {
  // Write center densities for each species of each component
  for (size_t i = 0; i < sim->component_list.size(); i++) {
    Component *comp = sim->component_list[i];
    for (auto it = comp->rho_center_map.begin();
         it != comp->rho_center_map.end(); it++) {
      Component::Species_Type species = it->first;
      ArrayXd rho_center = it->second;
      char species_char = Component::species_enum_to_char(species);
      std::ostringstream file_name_ss;
      file_name_ss << "rho_" << comp->name << "_" << species_char;
      if (NPROCS > 1) {
        // Make the file name look like homopolymer_a_a.04.dat
        file_name_ss << "." << std::setw(2) << std::setfill('0') << RANK;
      }
      file_name_ss << ".dat";
      std::string file_name(file_name_ss.str());
      fs::path file_path = output_dir / file_name;
      write_one_file(file_path, rho_center);
    }
  }

  if (sim->iter == 0) {
    // Write convolution (smearing) function for each species
    for (auto it = sim->conv_function_map.begin();
         it != sim->conv_function_map.end(); it++) {
      Component::Species_Type species = it->first;
      ArrayXd conv_function = it->second;
      char species_char = Component::species_enum_to_char(species);
      std::ostringstream file_name_ss;
      file_name_ss << "conv_func_" << species_char;
      if (NPROCS > 1) {
        // Make the file name look like conv_func_homopolymer_a_a.04.dat
        file_name_ss << "." << std::setw(2) << std::setfill('0') << RANK;
      }
      file_name_ss << ".dat";
      std::string file_name(file_name_ss.str());
      fs::path file_path = output_dir / file_name;
      write_one_file(file_path, conv_function);
    }

    // Write pair potentials
    for (int i = 0; i < Component::max_n_species; i++) {
      for (int j = i; j < Component::max_n_species; j++) {
        if (sim->pair_potential_arrays[i][j].size() == 0) {
          continue;
        }
        std::ostringstream file_name_ss;
        file_name_ss << "pair_potential_" << char('a' + i) << "_"
                     << char('a' + j);
        if (NPROCS > 1) {
          // Make the file name look like pair_potential_a_a.04.dat
          file_name_ss << "." << std::setw(2) << std::setfill('0') << RANK;
        }
        file_name_ss << ".dat";
        std::string file_name(file_name_ss.str());
        fs::path file_path = output_dir / file_name;
        write_one_file(file_path, sim->pair_potential_arrays[i][j]);
      }
    }
  }
}