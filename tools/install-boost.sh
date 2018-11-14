#!/bin/bash
#
# install-boost.sh
#
# Copyright (c) 2018 Ben Lindsay <benjlindsay@gmail.com>

proj_root="$(pwd)"
major="1"
minor="67"
patch="0"
version="$major.$minor.$patch"
boost_dir="boost_"$major"_"$minor"_"$patch
tarball="$boost_dir.tar.bz2"
url="https://dl.bintray.com/boostorg/release/$version/source/$tarball"

if [ ! -f /tmp/$tarball ]; then
  echo "downloading $url"
  if [ $(command -v wget) ]; then
    wget -O /tmp/$tarball $url
  elif [ $(command -v curl) ]; then
    curl -o /tmp/$tarball $url
  fi
else
  echo "/tmp/$tarball already downloaded."
fi

echo "Removing /tmp/$boost_dir..."
rm -rf /tmp/$boost_dir
cd /tmp
echo "Extracting /tmp/$tarball..."
tar -xjf $tarball
cd $boost_dir
echo "Running './bootstrap.sh'..."
./bootstrap.sh --prefix=$proj_root --with-libraries=filesystem
echo "Running './b2 install'..."
b2_output=$(./b2 install)
echo "$b2_output" | head -50
echo "...\nSkipping a ton of lines of './b2 install' output...\n..."
echo "$b2_output" | tail -50
echo "Removing *.so* and *.dylib files..."
rm -f $proj_root/lib/libboost*.so* $proj_root/lib/libboost*.dylib
