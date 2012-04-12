// -*- C++ -*- $Id$

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
