// summary_output.hpp
//
// Copyright (c) 2018 Ben Lindsay <benjlindsay@gmail.com>

#ifndef SUMMARY_OUTPUT_HPP
#define SUMMARY_OUTPUT_HPP

#include "../output.hpp"

class Summary_Output : public Output {
 public:
  Summary_Output(Sim *sim, std::vector<std::string> var_list, int print_freq,
                 fs::path file_path, int column_width, bool write_header);
  Summary_Output(Sim *sim, fs::path output_dir,
                 std::vector<std::string> options_list);
  virtual ~Summary_Output();
  virtual bool is_time_to_write();
  virtual void write();
  static const std::string default_file_name;

 private:
  void init(std::vector<std::string> options_list, int print_freq,
            fs::path file_path, int column_width, bool write_header);
  void init(fs::path output_dir, std::vector<std::string> options_list);
  int print_freq;
  int column_width;
  int iter_column_width;
  bool write_header;
  std::vector<std::string> options_list;
  fs::path file_path;
  fs::ofstream file;
};

#endif  // SUMMARY_OUTPUT_HPP
