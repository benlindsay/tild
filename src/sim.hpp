// sim.hpp
//
// Copyright (c) 2017 Ben Lindsay <benjlindsay@gmail.com>

#ifndef SIM_HPP
#define SIM_HPP

#include <string>  // std::string
#include <vector>  // std::vector
#include "Eigen/Dense"
#include "yaml-cpp/yaml.h"
#ifdef MPI
#include "fftw3-mpi.h"
#include "mpi.h"
#else
#include "fftw3.h"
#endif

#include "component.hpp"
#include "component_factory.hpp"
#include "globals.hpp"
#include "grid_utils.hpp"
#include "output.hpp"
#include "output_factory.hpp"
#include "utils.hpp"

class Output;

using Eigen::ArrayXd;   // Dynamically sized double Array
using Eigen::ArrayXcd;  // Dynamically sized double Array
using Eigen::ArrayXi;   // Dynamically sized int Array
using Eigen::ArrayXXi;  // Dynamically sized 2D int Array

// class FFTW_Utils;

class Sim {
 public:
  Sim(YAML::Node input);
  virtual ~Sim();
  virtual std::string get_var_as_string(std::string var_name, int str_len) = 0;
  virtual void init_default_summary_var_list() = 0;
  virtual void init_output_list(YAML::Node input) = 0;
  virtual void calculate_grid_densities();
  virtual void calculate_forces() = 0;
  virtual void move_particles() = 0;
  virtual void run();
  int get_global_index(int ix_global, int iy_global);
  int get_global_index(int ix_global, int iy_global, int iz_global);
  void init_component_list(YAML::Node input);
  void write_outputs();
  void init_fftw();
  void fftw_fwd(ArrayXd &in_array, ArrayXcd &out_array);
  void fftw_back(ArrayXcd &in_array, ArrayXd &out_array);
  void convolve(ArrayXd &array_1, ArrayXd &array_2, ArrayXd &convolved);

  std::string description;

  // FFTW_Utils *fftw_utils;

  int dim;
  int iter;
  int max_iter;
  double rho_0;
  double bond_length;
  double monomer_size;

  // Grid/box variables
  ArrayXd Lx;
  ArrayXi Nx;
  ArrayXd dx;
  ArrayXXi local_grid_indices;
  ArrayXXd local_grid_coords;
  ArrayXXd local_grid_k_coords;
  ArrayXd local_grid_k_magnitude;
  double V;
  int M;
  double grid_point_volume;
  int ML;
  int mesh_order;
  int n_subgrid_points;
  ArrayXXi weight_subgrid_index_shifts;

  // Potential variables
  std::vector<std::vector<ArrayXd> > pair_potential_arrays;

  std::vector<Output*> output_list;
  std::vector<Component*> component_list;
  std::vector<std::string> default_summary_var_list;

  // FFTW variables
  fftw_plan forward_plan;
  fftw_plan backward_plan;
  fftw_complex *fftw_in_array;
  fftw_complex *fftw_out_array;

 private:
  void init_box_vars(YAML::Node input);
};

#endif  // SIM_HPP
