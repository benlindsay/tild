// component_factory.cpp
//
// Copyright (c) 2018 Ben Lindsay <benjlindsay@gmail.com>

#include "component_factory.hpp"
#include "component_types/homopolymer.hpp"

Component* Component_Factory::New_Component(Sim* sim, YAML::Node input,
                                            YAML::Node component_params) {
  std::string component_type = component_params["type"].as<std::string>();
  if (component_type == "homopolymer") {
    Component::Species_Type species = Component::A;
    double vol_frac = -1;
    int n_segments_per_molecule = -1;
    for (YAML::const_iterator it = component_params.begin();
         it != component_params.end(); ++it) {
      std::string key = it->first.as<std::string>();
      YAML::Node value = it->second;
      if (key == "type") {
        continue;
      } else if (key == "species") {
        char species_char = value.as<char>();
        species = Component::species_char_to_enum(species_char);
      } else if (key == "vol_frac") {
        vol_frac = value.as<double>();
      } else if (key == "n_segments_per_molecule") {
        n_segments_per_molecule = value.as<int>();
      } else {
        utils::die("Can't recognize summary output parameter '" +
                   value.as<std::string>() + "'");
      }
    }
    if (vol_frac < 0) {
      utils::die("vol_frac for " + component_type + " is undefined");
    }
    if (n_segments_per_molecule < 0) {
      utils::die("n_segments_per_molecule for " + component_type +
                 " is undefined");
    }
    return new Homopolymer(sim, vol_frac, n_segments_per_molecule, species);
  } else {
    utils::die("Can't recognize component type '" + component_type + "'!");
    return NULL;
  }
  // Suppress "control may reach end of non-void function" warning
  return NULL;
}