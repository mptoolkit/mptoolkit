// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/openmp.cpp
//
// Copyright (C) 2022 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "common/openmp.h"
#include "common/environment.h"
#include <iostream>

#if defined(MULTITHREAD)

namespace omp
{

int MpNumThreads = getenv_or_default("MP_NUM_THREADS", -1);

void initialize(int Verbose)
{
   if (Verbose > 0)
   {
      std::cout << "Using OpenMP\n";
   }
   int NumProcs = omp_get_num_procs();
   if (MpNumThreads > -1)
   {
      if (Verbose > 0)
      {
         std::cout << "MP_NUM_THREADS is " << MpNumThreads << '\n';
      }
      omp_set_num_threads(MpNumThreads);
   }
   if (Verbose > 0)
   {
      std::cout << "Number of processors: " << NumProcs << '\n';
      std::cout << "Max threads per section: " << omp_get_max_threads() << '\n';
      #pragma omp parallel
      {
         #pragma omp master
         std::cout << "Actual number of threads: " << omp_get_num_threads() << '\n';
      }
      //std::cout << "Max threads: " << omp_get_thread_limit() << '\n';
   }
}

int threads_to_use(int Request)
{
   if (MpNumThreads == -1)
   {
      MpNumThreads = omp_get_max_threads();
   }
   return std::min(MpNumThreads, Request);
}

} // namespace omp

#endif
