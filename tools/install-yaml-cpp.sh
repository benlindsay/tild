#!/bin/bash
#
# install-yaml-cpp.sh
#
# Copyright (c) 2018 Ben Lindsay <benjlindsay@gmail.com>

YAML_CPP=/tmp/yaml-cpp
PROJECT_ROOT=$(git rev-parse --show-toplevel)

rm -rf $YAML_CPP
git clone https://github.com/jbeder/yaml-cpp $YAML_CPP && \
  git checkout -b 2018-05-14 4fb1c4b92bf8d9 && \
  mkdir $YAML_CPP/build && \
  cd $YAML_CPP/build && \
  CC=$(which gcc) CXX=$(which g++) \
    cmake -DCMAKE_INSTALL_PREFIX=$PROJECT_ROOT .. && \
  make && \
  make install && \
  rm -rf $PROJECT_ROOT/include/gmock  $PROJECT_ROOT/include/gtest
