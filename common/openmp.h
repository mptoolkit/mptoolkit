// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// common/openmp.h
//
// Copyright (C) 2012-2021 Ian McCulloch <ian@qusim.net>
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

#if !defined(MPTOOLKIT_COMMON_OPENMP_H)
#define MPTOOLKIT_COMMON_OPENMP_H

#include "config.h"

#if defined(HAVE_OPENMP)
   #if !defined(MULTITHREAD)
      #error "Build configuration error: OpenMP support requires MULTITHREAD."
   #endif
   #include "openmp_thread.h"
#else
   #if defined(_OPENMP)
      #error "Build configuration error: _OPENMP is defined, but the build is configured without OpenMP!  Rerun the configure --with-openmp, or disable OpenMP compiler."
   #else
      #include "openmp_dummy.h"
   #endif
#endif

#endif
