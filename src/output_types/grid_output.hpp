// grid_output.hpp
//
// Copyright (c) 2018 Ben Lindsay <benjlindsay@gmail.com>

#ifndef GRID_OUTPUT_HPP
#define GRID_OUTPUT_HPP

#include "../output.hpp"

class Grid_Output : public Output {
 public:
  Grid_Output(Sim *sim, fs::path output_dir, int print_freq, int column_width,
              bool write_header, bool gnuplot_compatible);
  Grid_Output(Sim *sim, fs::path output_dir);
  virtual ~Grid_Output(){};
  virtual bool is_time_to_write();
  virtual void write();
  virtual void write_iter_0();

 private:
  void init(fs::path output_dir, int print_freq, int column_width,
            bool write_header, bool gnuplot_compatible);
  void init(fs::path output_dir);
  void write_one_file(fs::path file_path, ArrayXd& data);
  fs::path output_dir;
  int print_freq;
  int column_width;
  bool write_header;
  bool gnuplot_compatible;
  fs::ofstream file;
};

#endif  // GRID_OUTPUT_HPP