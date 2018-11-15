// lammpstrj_output.hpp
//
// Copyright (c) 2018 Ben Lindsay <benjlindsay@gmail.com>

#ifndef LAMMPSTRJ_OUTPUT_HPP
#define LAMMPSTRJ_OUTPUT_HPP

#include "../output.hpp"

class Lammpstrj_Output : public Output {
 public:
  Lammpstrj_Output(Sim *sim, std::vector<Component *> component_list,
                   fs::path output_dir, int print_freq, std::string name,
                   bool one_frame_per_file);
  virtual ~Lammpstrj_Output();
  virtual bool is_time_to_write();
  virtual void write_iter_0();
  virtual void write();
  static const std::string default_name;

 private:
  void init(std::vector<Component *> component_list, fs::path output_dir,
            int print_freq, std::string name, bool one_frame_per_file);
  void write_one_frame();
  std::vector<Component *> component_list;
  int print_freq;
  std::string name;
  bool one_frame_per_file;
  fs::path file_path;
  fs::ofstream file;
};

#endif  // LAMMPSTRJ_OUTPUT_HPP