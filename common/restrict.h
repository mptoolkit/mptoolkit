/* -*- C++ -*- $Id$
  restrict.h

  Makes the 'restrict' keyword available, if the compiler supports it.

  Created 2003-10-18 Ian McCulloch
*/

#if !defined(RESTRICT_H_DFJIROEUIOJOI)
#define RESTRICT_H_DFJIROEUIOJOI

#if defined(restrict)
// restrict already defined, do nothing

#elif defined(__KCC)
// restrict exists in KCC, no need to do anything special

#elif defined(__sgi)
// restrict in SGI is the double underscore version
#define restrict __restrict

#elif defined(__INTEL_COMPILER)
// restrict exists in icc, no need to do anything special

#elif defined(__GNUC__)
// double-underscore restrict always works in gcc
#define restrict __restrict__

#else
// no restrict :(
#define restrict
#endif

#endif
