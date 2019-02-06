// canonical_sim.hpp
//
// Copyright (c) 2017 Ben Lindsay <benjlindsay@gmail.com>

#ifndef CANONICAL_SIM_HPP
#define CANONICAL_SIM_HPP

#include "../component.hpp"
#include "../sim.hpp"

class Canonical_Sim : public Sim {
 public:
  virtual void init(YAML::Node input);
  virtual ~Canonical_Sim(){};
  virtual void init_default_summary_options_list();
  virtual void init_output_list(YAML::Node input);
};

#endif  // CANONICAL_SIM_HPP
