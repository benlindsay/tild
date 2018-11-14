#!/bin/bash
#
# install-clang-format.sh
#
# Copyright (c) 2018 Ben Lindsay <benjlindsay@gmail.com>

if [ $# -gt 0 ]; then
  REPLY="$1"
fi

PROJECT_ROOT=$(git rev-parse --show-toplevel)
CLANG_FORMAT="$PROJECT_ROOT/tools/clang-format"

if [ -f $CLANG_FORMAT ]; then
  if [ $($CLANG_FORMAT --version | grep '3.8.0' | wc -l) -ge 1 ]; then
    echo "You already have tools/clang-format version 3.8.0."
    exit 0
  else
    echo "You have tools/clang-format, but it isn't version 3.8.0."
  fi
fi

echo 'Here are the results of $(which clang-format):'
echo $(which clang-format)

if [ $(clang-format --version | grep '3.8.0' | wc -l) -ge 1 ]; then
  cf_path=$(which clang-format)
  cf_dir=$(dirname $cf_path)
  echo "You already have clang-format 3.8.0 installed in $cf_dir."
  echo "Copying to tools/ :"
  cp $cf_path $PROJECT_ROOT/tools/clang-format
  exit 0
fi

# These options come from http://releases.llvm.org/download.html, v 3.8.0
# Any option after that won't let you just grab a pre-built clang-format
# executable.
options=(
  "Mac OS X"
  "Fedora23 x86_64 Linux"
  "CentOS 6 x86_64"
  "x86_64 Ubuntu 14.04"
)
urls=(
  "http://releases.llvm.org/3.8.0/clang+llvm-3.8.0-x86_64-apple-darwin.tar.xz"
  "http://releases.llvm.org/3.8.0/clang+llvm-3.8.0-x86_64-fedora23.tar.xz"
  "http://releases.llvm.org/3.8.0/clang+llvm-3.8.0-linux-x86_64-centos6.tar.xz"
  "http://releases.llvm.org/3.8.0/clang+llvm-3.8.0-x86_64-linux-gnu-ubuntu-14.04.tar.xz"
)

echo "I can't find clang-format version 3.8.0. Let's install it."

if [ -z $REPLY ]; then

  title="Which OS is closest to yours?"
  prompt="Pick an option:"

  echo "$title"
  PS3="$prompt "
  select opt in "${options[@]}" "Cancel"; do

    case "$REPLY" in

      1 ) break ;;
      2 ) break ;;
      3 ) break ;;
      4 ) break ;;
      $(( ${#options[@]}+1 )) ) echo "Goodbye!"; break ;;
      *) echo "Invalid option. Try another one."; continue ;;

    esac

  done
else
  opt=${options[$(($REPLY-1))]}
fi

echo "You picked $opt."
url=${urls[$(($REPLY-1))]}
tarball=$(basename $url)

if [ ! -f /tmp/$tarball ]; then
  echo "downloading $url"
  if [ $(command -v wget) ]; then
    wget -O /tmp/$tarball $url
  elif [ $(command -v curl) ]; then
    curl -o /tmp/$tarball $url
  fi
else
  echo "Already downloaded $tarball."
fi

src_dir=$(basename $tarball .tar.xz)
if [ ! -d /tmp/$src_dir ]; then
  tar xfJ /tmp/$tarball -C /tmp/
else
  echo "Already extracted $tarball."
fi

echo "Copying clang-format to tools/clang-format"
cp /tmp/$src_dir/bin/clang-format $PROJECT_ROOT/tools/clang-format

echo "Checking to make sure it's right:"
if [ $($CLANG_FORMAT --version | grep '3.8.0' | wc -l) -ge 1 ]; then
  echo "You've got version 3.8.0!"
  echo "I recommend copying clang-format to somewhere on your path if it's not"
  echo "there already."
  echo "Try something like 'cp tools/clang-format ~/bin/'"
  exit 0
else
  echo "Uh-oh, something went wrong and you don't have version 3.8.0..."
  exit 1
fi
