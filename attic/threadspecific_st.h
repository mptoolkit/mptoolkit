// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// attic/threadspecific_st.h
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
  single-threaded (dummy) version of threadspecific.h

  Created 2002-10-26 Ian McCulloch
*/

#if defined(MULTITHREAD)
#error "Cannot use single thread header threadspecific_st.h with MULTITHREAD."
#endif

namespace pthread
{

template <class T>
class thread_specific
{
   public:
      thread_specific() {}
      thread_specific(T const& Init) : data(Init) {}

      operator T&() { return data; }
      operator T const&() const { return data; }

   private:
      thread_specific(thread_specific<T> const&); // not implemented
      thread_specific<T>& operator=(thread_specific<T> const&); // not implemented

      T data;
};

} // namespace pthread
