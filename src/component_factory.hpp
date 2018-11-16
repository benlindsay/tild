// component_factory.hpp
//
// Copyright (c) 2018 Ben Lindsay <benjlindsay@gmail.com>

#ifndef COMPONENT_FACTORY_HPP
#define COMPONENT_FACTORY_HPP

#include "component.hpp"
#include "component_types/homopolymer.hpp"
#include "yaml-cpp/yaml.h"

class Component_Factory {
 public:
  static Component *New_Component(Sim *sim, YAML::Node input,
                                  YAML::Node component_params);
};

#endif  // COMPONENT_FACTORY_HPP
