#!/bin/bash
#
# install-yaml-cpp.sh
#
# Copyright (c) 2018 Ben Lindsay <benjlindsay@gmail.com>

YAML_CPP=/tmp/yaml-cpp
PROJECT_ROOT=$(git rev-parse --show-toplevel)

rm -rf $YAML_CPP
git clone https://github.com/jbeder/yaml-cpp $YAML_CPP && \
  cd $YAML_CPP && \
  git checkout -b 2018-05-14 4fb1c4b92bf8d94b32ebccdd890407d45b3bc794 && \
  mkdir build && \
  cd build && \
  CC=$(which gcc) CXX=$(which g++) cmake -DCMAKE_INSTALL_PREFIX=$PROJECT_ROOT \
    -DYAML_CPP_BUILD_TESTS=OFF .. && \
  make && \
  make install
