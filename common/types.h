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

// type aliases for converting between real and complex numbers.
// Provides:
// real_t<T>    : real-valued type associated with T
// complex_t<T> : complex-valued type associated with T
// scalar_t<T>  : scalar-valued type associated with T
//
// if T is already a scalar, then scalar_t<T> is the identity.
// To use: specialize the template ScalarTypes<T> to and define the members
// real_t, complex_t and scalar_t.

template <typename T>
struct ScalarTypes
{
   using real_t = T;
   using complex_t = std::complex<T>;
   using scalar_t = T;
};

template <typename T>
struct ScalarTypes<std::complex<T>> : ScalarTypes<T> {};

template <typename T>
using real_t = typename ScalarTypes<T>::real_t;

template <typename T>
using complex_t = typename ScalarTypes<T>::complex_t;

template <typename T>
using scalar_t = typename ScalarTypes<T>::scalar_t;

#if defined(HAVE_FLOAT128)
#include <boost/multiprecision/float128.hpp>
using float128 = boost::multiprecision::float128;
using complex256 = std::complex<float128>;
#endif

using float64 = double;
using complex128 = std::complex<double>;

using float32 = float;
using complex64 = std::complex<float>;

#if defined(USE_FLOAT128)
using real = float128;
using complex = complex256;
#else
using real = float64;
using complex = complex128;
#endif

#endif
