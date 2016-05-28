// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/mutex.h
//
// Copyright (C) 2002-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
  mutex.h

  Split the locking primitives out of thread.h, to allow for
  a dummy non-threaded implementation.

  Created 2002-07-22 Ian McCulloch
*/

#if !defined(MUTEX_H_HFUIY5T78Y78FHUIH789T7F43YF89ER8YH5RAHY)
#define MUTEX_H_HFUIY5T78Y78FHUIH789T7F43YF89ER8YH5RAHY

#if defined(MULTITHREAD)
#include "mutex_mt.h"
#else
#include "mutex_st.h"
#endif

#endif
