// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/fortran.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Reseach publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

/*
  fortran.h

  Interface macros and types for calling fortran from C++.
  The F77NAME is somewhat platform dependent.

  General guidelines for calling fortran from C++:
  The typical procedure we use is a 2-level scheme, where the inner namespace 'raw' contains
  the original extern "C" prototype, and then an outer level function of the same
  name does some basic fixes to the types etc.  In this outer scope function,
  all non-array parameters should be pass by value rather than pass by pointer.
  Input arrays that are not modified by fortran should be passed by const pointer (double const* Ptr),
  (should these also be restrict'ed?)
  Output arrays that are modified by fortran should be passed by restrict pointer (double* restrict Ptr).
  
  Handling of complex types is a bit tricky, for a number of reasons.  C++98 contains no builtin
  complex type, but on any platform where we have a hope of interfacing C++ with fortran directly,
  it is likely that the C++ complex<double> type will be layout compatible with
  fortran COMPLEX*16.  Thus it should always be possible to reinterpret_cast to/from
  fortran complex types.  There is a Fortran::complex type defined here, but I don't recommend using it.

  Fortran functions that return a value of type complex are tricky to handle as the
  calling convention for this is not fixed.  On linux/x86, a pointer to the return value
  is passed as the first parameter to the function.  On amd64, the platform ABI specifies 
  that complex return values are passed in registers, so no special handling is needed.  
  To handle these two cases, define one (and only one) of the symbols
  FORTRAN_COMPLEX_RETURN_FIRST_ARG or FORTRAN_COMPLEX_RETURN_IN_REGISTER.
  If any of these symbols are defined, then this header also defines
  HAVE_FORTRAN_COMPLEX_RETURN.  Wrapper functions can make use of these symbols to
  declare the correct prototypes.  See blas1f.h for an example: the prototype is only declared
  if HAVE_FORTRAN_COMPLEX_RETURN is defined (otherwise the compex-returning functions
  are not available), and in the body of the wrapper function, the calling convention is tested
  to see how to make the forward call to the fortran function.
*/

#if !defined(FORTRAN_H_SFDHH348U9J3JCERJEFJIFUCFUJFU803FCJ)
#define FORTRAN_H_SFDHH348U9J3JCERJEFJIFUCFUJFU803FCJ

#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#include <complex>

#if defined(FORTRAN_DOUBLE_UNDERSCORE)
#define F77NAME(x) x##__
#else
#define F77NAME(x) x##_
#endif

#if defined(FORTRAN_COMPLEX_RETURN_FIRST_ARG) || defined(FORTRAN_COMPLEX_RETURN_IN_REGISTER)
#if !defined(HAVE_FORTRAN_COMPLEX_RETURN)
#define HAVE_FORTRAN_COMPLEX_RETURN
#endif
#endif

// fortran typedefs

namespace Fortran
{

typedef int integer;
typedef int logical;
typedef double real;

struct complex
{
   double real;
   double imag;

   operator std::complex<double>() const { return std::complex<double>(real, imag); }
};

} // namespace

#endif
