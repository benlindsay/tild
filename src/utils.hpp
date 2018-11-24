// utils.hpp
//
// Copyright (c) 2017 Ben Lindsay <benjlindsay@gmail.com>

#ifndef UTILS_HPP
#define UTILS_HPP

#include <algorithm>  // std::transform
#include <fstream>    // std::ofstream
#include <iostream>   // std::cout, std::endl
#include <string>     // std::string
#include "Eigen/Dense"
#include "globals.hpp"
#include "yaml-cpp/yaml.h"

#define MAX_PATH_LEN 256

using Eigen::ArrayXXd;  // Dynamically sized 2D double Array
using Eigen::ArrayXd;   // Dynamically sized 2D double Array

namespace utils {

void die();
void die(std::string);
void die(std::stringstream&);
bool check_pbc(ArrayXXd& array, ArrayXd& Lx);
void enforce_pbc(ArrayXXd& array, ArrayXd& Lx);
void enforce_pbc(ArrayXXd& array, ArrayXd& Lx, int i_site);
void mpi_init_wrapper(int, const char**);
void mpi_finalize_wrapper();
void print_one_line(std::string);
void print_one_line(const char*);
void print_one_line(std::stringstream&);
void print_one_line(std::ofstream&, std::string);
void print_one_line(std::ofstream&, std::stringstream&);
std::string to_lower(std::string);
void to_lower(YAML::Node);

}  // namespace utils

#endif  // UTILS_HPP
