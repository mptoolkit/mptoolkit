// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/openmp.h
//
// Copyright (C) 2017 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

//
// Dummy implementation of OpenMP runtime functions
//

#if !defined(MPTOOLKIT_COMMON_OPENMP_DUMMY_H)
#define MPTOOLKIT_COMMON_OPENMP_DUMMY_H

namespace omp
{

inline
void initialize()
{
}

inline
int threads_to_use(int Request)
{
   return 1;
}

template <typename T>
T
parallel_sum(std::vector<T>&& x)
{
   CHECK(x.size() > 0);

   T Result = x[0];
   for (unsigned n = 1; n < x.size(); ++n)
   {
      Result += x[n];
   }
   return Result;
}

} // namespace omp

#endif
