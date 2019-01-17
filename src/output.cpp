// output.cpp
//
// Copyright (c) 2019 Ben Lindsay <benjlindsay@gmail.com>

#include "output.hpp"

void Output::write_double(std::stringstream &ss, double value,
                          int string_length) {
  // Assuming we want 2 characters go to "#." and 4 characters go to exponent,
  // i.e. "e+00", string length - 6 gives the number
  // of digits we can put after the decimal
  int scientific_precision = string_length - 7;
  // Set number of digits after decimal point in double notation
  ss.precision(scientific_precision);
  // Set minimum width of output string (left-padded with spaces)
  ss << std::setw(string_length);
  // Add scientific notation flag to force specified number of digits after
  // decimal point
  ss << std::scientific;
  ss << value;
}
