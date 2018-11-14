#!/bin/bash
#
# clang-format-diff.sh
#
# Copyright (c) 2018 Ben Lindsay <benjlindsay@gmail.com>

PROJECT_ROOT=$(git rev-parse --show-toplevel)
CLANG_FORMAT="$PROJECT_ROOT/tools/clang-format"

if [ ! -f $CLANG_FORMAT ]; then
  echo "$CLANG_FORMAT is not a valid file. Run 'make clang-format' then retry."
  exit 1
fi

cd $PROJECT_ROOT

file_list=$( find src test ! -path 'test/catch.hpp' \
  \( -name '*.cpp' -or -name '*.hpp' -or -name '*.c' -or -name '*.h' \) \
)

# Show the changes made by diffing all source files except catch.hpp.
# Failing exit code if any changes are made.

diff -u <(cat $file_list) <($CLANG_FORMAT -style=file $file_list)
