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
#include <string>

namespace omp
{

inline
void initialize(int Verbose = 0)
{
   int NumProcs = omp_get_num_procs();
   char const* Str = getenv("MP_NUM_THREADS");
   if (Str)
   {
      int n = std::stoi(std::string(Str));
      omp_set_num_threads(n);
      if (Verbose > 0)
      {
         std::cout << "MP_NUM_THREADS is " << n << '\n';
      }
   }
   if (Verbose > 0)
   {
      std::cout << "Number of processors: " << NumProcs << '\n';
      #pragma omp parallel
      {
         #pragma omp single
         std::cout << "Number of threads: " << omp_get_num_threads() << '\n';
      }
      std::cout << "Max threads per section: " << omp_get_max_threads() << '\n';
      std::cout << "Max threads: " << omp_get_thread_limit() << '\n';
   }
}

inline
int threads_to_use(int Request)
{
   return Request;
}

} // namespace omp

#endif
