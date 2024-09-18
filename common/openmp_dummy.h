// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// common/openmp_dummy.h
//
// Copyright (C) 2015-2022 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2022 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

//
// Dummy implementation of OpenMP runtime functions
//

#if !defined(MPTOOLKIT_COMMON_OPENMP_DUMMY_H)
#define MPTOOLKIT_COMMON_OPENMP_DUMMY_H

#include <iostream>

namespace omp
{

inline
void initialize(int Verbose = 0)
{
   if (Verbose > 0)
   {
      std::cout << "OpenMP is disabled.\n";
   }
}

inline
int threads_to_use(int Request)
{
   return 1;
}

} // namespace omp

#endif
