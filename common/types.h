// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// common/types.h
//
// Copyright (C) 2017 Ian McCulloch <ian@qusim.net>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

/*
  floating point typedefs.

  Defines types in the global namespace:

  float64       64-bit floating point type
  complex128    complex number made from two 64-bit numbers for the real and imaginary parts

  We also define typedefs for the default floating point type.

  real == float64
  complex == complex128
*/

#if !defined(MPTOOLKIT_COMMON_TYPES_H)
#define MPTOOLKIT_COMMON_TYPES_H

#include <complex>

#if defined(HAVE_FLOAT128) || defined(USE_FLOAT128)
#error "Boost.Multiprecision float128 support has been removed from MPTK."
#endif

using real64 = double;
using complex128 = std::complex<double>;

using real = real64;
using complex = complex128;

#endif
