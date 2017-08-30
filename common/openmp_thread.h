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

#if !defined(MPTOOLKIT_COMMON_OPENMP_THREAD_H)
#define MPTOOLKIT_COMMON_OPENMP_THREAD_H

#include <omp.h>

namespace omp
{

inline
void initialize()
{
   omp_set_dynamic(true);
   omp_set_nested(true);
   int n = omp_get_thread_limit();
   omp_set_num_threads(n); // seems we need to initialize this???
#pragma omp parallel num_threads(100)
   {
#pragma omp single
      std::cout << "Number of threads: " << omp_get_num_threads() << '\n';
   }
   std::cout << "Max threads per section: " << omp_get_max_threads() << '\n';
   std::cout << "Max threads: " << omp_get_thread_limit() << '\n';
}

inline
int threads_to_use(int Request)
{
   return Request;
}

// some parallel algorithms
template <typename T>
T
parallel_sum(std::vector<T>&& x)
{
   CHECK(x.size() > 0);
#if 0
   T Result = x[0];
   for (unsigned n = 1; n < x.size(); ++n)
   {
      Result += x[n];
   }
   return Result;

#else
   int TotalSz = x.size();
   int Stride = 1;
   while (Stride < TotalSz)
   {
      int Sz = TotalSz / (2*Stride);
#pragma omp parallel for schedule(dynamic) num_threads(omp::threads_to_use(Sz))
      for (int n = 0; n < Sz; ++n)
      {
         x[n*Stride*2] += x[n*Stride*2 + Stride];
      }
      Stride *= 2;
   }
   return x[0];
#endif
}

} // namespace omp

#endif
