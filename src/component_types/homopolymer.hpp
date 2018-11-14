// homopolymer.hpp
//
// Copyright (c) 2018 Ben Lindsay <benjlindsay@gmail.com>

#ifndef HOMOPOLYMER_HPP
#define HOMOPOLYMER_HPP

#include "../component.hpp"

class Homopolymer : public Component {
 public:
  Homopolymer(Sim *sim, int n_segments, Component::Species_Type species);
  virtual ~Homopolymer(){};
  int n_segments;
};

#endif  // HOMOPOLYMER_HPP