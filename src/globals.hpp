// globals.hpp
//
// Copyright (c) 2017 Ben Lindsay <benjlindsay@gmail.com>

#ifndef GLOBALS_HPP
#define GLOBALS_HPP

// #include <gperftools/profiler.h>
#include <boost/filesystem.hpp>  // boost::filesystem::path,canonical,absolute
#include <complex>               // std::complex

namespace fs = boost::filesystem;

#ifdef MAIN_HPP
int RANK = 0;
#else
extern int RANK;
#endif

#ifdef MAIN_HPP
int NPROCS = 1;
#else
extern int NPROCS;
#endif

#ifdef MAIN_HPP
std::complex<double> I(0.0, 1.0);
#else
extern std::complex<double> I;
#endif

#define PI 3.14159265358979323846263383

#endif  // GLOBALS_HPP
