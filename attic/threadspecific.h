// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// attic/threadspecific.h
//
// Copyright (C) 2012-2016 Ian McCulloch <ian@qusim.net>
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
  threadspecific.h

  Created 2002-10-26 Ian McCulloch

  Handles thread-specific data.  Where possible, it is better to store thread-specific
  data as members in the class derived from thread<Arg, Result>, however
  thread-specific data is sometimes necessary in library code that
  doesn't have access to the thread objects.
  To make thread-specific data of type T, create a single
  object of type thread_specific<T> which is visible to all threads.
  Then each thread that calls thread_specific<T>::value() will
  get a distinct object of type T.  The thread-specific
  object is default constructed on the first call to thread_specific<T>::value()
  by each thread.  The object is properly destructed then a thread is destroyed.
  Note the pthread limits on the number of thread-specific data entries.
*/

#if !defined(THREADSPECIFIC_H_SDHCE4UYGU53YT87HJCIUO24H897HPAS)
#define THREADSPECIFIC_H_SDHCE4UYGU53YT87HJCIUO24H897HPAS

#if defined(MULTITHREAD)
#include "threadspecific_mt.h"
#else
#include "threadspecific_st.h"
#endif

#endif
