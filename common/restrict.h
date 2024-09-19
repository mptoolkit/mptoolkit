// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// common/restrict.h
//
// Copyright (C) 2003-2016 Ian McCulloch <ian@qusim.net>
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
