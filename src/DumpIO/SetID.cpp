#include <cassert>
#include <iostream>
#include <vector>
#include "DumpIO.hpp"

void DumpIO::SetID(std::vector<int>& id_t_) {
  if (run_input or !run_output) {
    std::cout << "ERROR: running input function after opening output";
    assert(!run_input or run_output);
  }

  id_t = id_t_;
  set_id = true;

  if (set_pos) {
    if (size_t(n_p) != id_t_.size()) {
      std::cout << "Number of particles mismatch!\n";
      assert(size_t(n_p) == id_t_.size());
    }
  } else {
    n_p = id_t_.size();
  }
}
