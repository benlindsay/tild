// semi_grand_sim.hpp
//
// Copyright (c) 2019 Ben Lindsay <benjlindsay@gmail.com>

#ifndef SEMI_GRAND_SIM_HPP
#define SEMI_GRAND_SIM_HPP

#include <vector>
#include "../component.hpp"
#include "../sim.hpp"
#include "../utils.hpp"

class Semi_Grand_Sim : public Sim {
 public:
  virtual void init(YAML::Node input) override;
  virtual ~Semi_Grand_Sim(){};
  virtual void init_default_summary_options_list();
  virtual void init_output_list(YAML::Node input);
  virtual void init_component_list(YAML::Node input) override;
  virtual void run() override;
  void update_fractional_presence();
  double calculate_dU_chi_kappa_dlambda_1_m_1(Component *comp_1,
                                              Component *comp_2);

  double partial_step_rate;
  std::vector<int> swappable_components_list;
};

#endif  // SEMI_GRAND_SIM_HPP
