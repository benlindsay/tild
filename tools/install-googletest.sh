#!/usr/bin/env bash
#
# install-googletest.sh
#
# Copyright (c) 2019 Ben Lindsay <benjlindsay@gmail.com>

GTEST_REPO=/tmp/googletest
GTEST_DIR=$GTEST_REPO/googletest
GMOCK_DIR=$GTEST_REPO/googlemock
PROJECT_ROOT=$(git rev-parse --show-toplevel)
TEST_BUILDDIR=testbuild

rm -rf $GTEST_REPO
git clone https://github.com/google/googletest $GTEST_REPO && \
  cd $GTEST_REPO && \
  git checkout -b 2019-01-14 eb9225ce361affe561592e0912320b9db84985d0 && \
  rm -rf $PROJECT_ROOT/include/gtest $PROJECT_ROOT/include/gmock && \
  cp -r $GTEST_DIR/include/gtest $PROJECT_ROOT/include/ && \
  cp -r $GMOCK_DIR/include/gmock $PROJECT_ROOT/include/ && \
  cd $PROJECT_ROOT && \
  mkdir -p $TEST_BUILDDIR && \
  g++ -std=c++11 -I $GTEST_DIR -I $GMOCK_DIR -I include -pthread \
    -c $GTEST_DIR/src/gtest-all.cc  -o $TEST_BUILDDIR/gtest-all.o && \
  g++ -std=c++11 -I $GTEST_DIR -I $GMOCK_DIR -I include -pthread \
    -c $GTEST_DIR/src/gtest_main.cc -o $TEST_BUILDDIR/gtest_main.o && \
  g++ -std=c++11 -I $GTEST_DIR -I $GMOCK_DIR -I include -pthread \
    -c $GMOCK_DIR/src/gmock-all.cc -o $TEST_BUILDDIR/gmock-all.o && \
  ar -rv lib/libgtest.a $TEST_BUILDDIR/gtest-all.o && \
  ar -rv lib/libgtest_main.a $TEST_BUILDDIR/gtest-all.o $TEST_BUILDDIR/gtest_main.o && \
  ar -rv lib/libgmock.a $TEST_BUILDDIR/gtest-all.o $TEST_BUILDDIR/gmock-all.o
