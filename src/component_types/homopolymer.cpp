// ft_homopolymer.cpp
//
// Copyright (c) 2018 Ben Lindsay <benjlindsay@gmail.com>

#include "homopolymer.hpp"

Homopolymer::Homopolymer(Sim *sim, int _n_segments,
                         Component::Species_Type species)
    : Component(sim) {
  n_segments = _n_segments;
}