// ENDHEADER

/*
  floating point typedefs.

  Defines types in the global namespace:

  float64       64-bit floating point type
  complex128    complex number made from two 64-bit numbers for the real and imaginary parts

  if HAVE_FLOAT128 is defined, then we use the boost multiprecision float128 type and define
  typedefs for

  float128      128-bit floating point type
  complex256    complex number made from two 128-bit numbers for the real and imaginary parts

  We also define some typedefs for a 'default' floating point type.
  If USE_FLOAT128 is defined, then

  real == float128
  complex == complex256

  otherwise we have

  real == float64
  complex == complex128
*/

#if !defined(MPTOOLKIT_COMMON_TYPES_H)
#define MPTOOLKIT_COMMON_TYPES_H

#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#include <complex>

#if defined(HAVE_FLOAT128)
#include <boost/multiprecision/float128.hpp>
using float128 = boost::multiprecision::float128;
using complex256 = std::complex<float128>;
#endif

using real64 = double;
using complex128 = std::complex<double>;

#if defined(USE_FLOAT128)
using real = float128;
using complex = complex256;
#else
using real = double;
using complex = std::complex<double>;
#endif

#endif
