// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// attic/threadspecific_mt.h
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

#if !defined(MULTITHREAD)
#error "Multithread header threadspecific_mt.h cannot be included without MULTITHREAD."
#endif

#include "threads.h"

namespace pthread
{

template <class T>
class thread_specific
{
   public:
      thread_specific();
      explicit thread_specific(T const& InitialValue);
      ~thread_specific();

      operator T&();
      operator T const&() const;

   private:
      thread_specific(thread_specific<T> const&); // not implemented
      thread_specific<T>& operator=(thread_specific<T> const&); // not implemented

      ::pthread_key_t key;
      T Init;
};

} // namespace pthread

#include "threadspecific_mt.cc"
