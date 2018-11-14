// utils.cpp
//
// Copyright (c) 2017 Ben Lindsay <benjlindsay@gmail.com>

#include "utils.hpp"

void utils::die() {
  std::string message = "";
  utils::die(message);
}

void utils::die(std::string message) {
  utils::print_one_line(message);
  utils::print_one_line("Exiting program.");

#ifdef MPI
  MPI_Finalize();
#endif

  exit(1);
}

void utils::die(std::stringstream& message) { utils::die(message.str()); }

void utils::mpi_init_wrapper(int argc, const char* argv[]) {
#ifdef MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &RANK);
  MPI_Comm_size(MPI_COMM_WORLD, &NPROCS);
#endif
}

void utils::mpi_finalize_wrapper() {
#ifdef MPI
  MPI_Finalize();
#endif
}

void utils::print_one_line(std::string line) {
  if (RANK == 0) {
    std::cout << line << std::endl;
  }
}

void utils::print_one_line(const char* line) {
  std::string line_str = line;
  utils::print_one_line(line_str);
}

void utils::print_one_line(std::stringstream& line) {
  utils::print_one_line(line.str());
}

void utils::print_one_line(std::ofstream& file, std::string line) {
  if (RANK == 0) {
    file << line << std::endl;
  }
}

void utils::print_one_line(std::ofstream& file, std::stringstream& line) {
  utils::print_one_line(file, line.str());
}

std::string utils::to_lower(std::string str) {
  // transform is from algorithm library
  // See link below:
  // stackoverflow.com/questions/313970/how-to-convert-stdstring-to-lower-case
  std::transform(str.begin(), str.end(), str.begin(), ::tolower);
  return str;
}

void utils::to_lower(YAML::Node node) {
  if (node.IsNull()) {
    return;
  }
  if (node.IsScalar()) {
    node = to_lower(node.as<std::string>());
  } else if (node.IsMap()) {
    for (YAML::iterator it = node.begin(); it != node.end(); ++it) {
      to_lower(it->first);
      to_lower(it->second);
    }
  } else if (node.IsSequence()) {
    for (std::size_t i = 0; i < node.size(); i++) {
      to_lower(node[i]);
    }
  }
}
