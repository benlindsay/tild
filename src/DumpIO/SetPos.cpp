#include <cassert>
#include <iostream>
#include "DumpIO.hpp"

void DumpIO::SetPos(std::vector<std::vector<double> >& pos_t_) {
  if (run_input or !run_output) {
    std::cout << "ERROR: running input function after opening output";
    assert(!run_input or run_output);
  }

  pos_t = pos_t_;
  set_pos = true;

  if (set_id) {
    if (size_t(n_p) != pos_t_.size()) {
      std::cout << "Number of particles mismatch!";
      assert(size_t(n_p) == pos_t_.size());
    }
  } else {
    n_p = pos_t_.size();
  }

  if (set_bb) {
    if (size_t(d) != pos_t_[0].size()) {
      std::cout << "Dimensions mismatch!";
      assert(size_t(d) == pos_t_[0].size());
    }
  } else {
    d = pos_t_[0].size();
  }
}
