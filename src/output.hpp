// output.hpp
//
// Copyright (c) 2018 Ben Lindsay <benjlindsay@gmail.com>

#ifndef OUTPUT_HPP
#define OUTPUT_HPP

#include <boost/filesystem/fstream.hpp>  // fs::ofstream
#include <iomanip>                       // std::setw, std::setprecision
#include <iostream>                      // std::cout, std::endl
#include "utils.hpp"

class Sim;

class Output {
 public:
  Output(Sim *sim) : sim(sim){};
  virtual ~Output(){};
  static void write_double(std::stringstream &ss, double value,
                           int string_length);
  virtual bool is_time_to_write() = 0;
  virtual void write() = 0;
  static const int default_column_width = 15;
  static const int default_print_freq = 50;

  Sim *sim;
};

#include "sim.hpp"

#endif  // OUTPUT_HPP
