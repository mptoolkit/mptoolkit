// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// quantumnumbers/test/testcoupling.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "quantumnumbers/coupling.h"
#include "common/trace.h"

int main()
{
   // check that the caching version of Racah gives the same result as the non-caching version.
   // This is therefore a test of the canonicalization of the 6j coefficients using the symmetry relations
   // and the hashing mechanism.
   half_int Max = 6;

   for (half_int j1 = 0; j1 <= Max; j1 += 0.5)
   {
      for (half_int j2 = 0; j2 <= Max; j2 += 0.5)
      {
	 for (half_int j3 = 0; j3 <= Max; j3 += 0.5)
	 {
	    for (half_int j4 = 0; j4 <= Max; j4 += 0.5)
	    {
	       for (half_int j5 = 0; j5 <= Max; j5 += 0.5)
	       {
		  if (!is_triangle(j1,j2,j5) || !is_triangle(j3,j4,j5))
		     continue;

		  for (half_int j6 = 0; j6 <= Max; j6 += 0.5)
		  {
		     if (!is_triangle(j1,j3,j6) || !is_triangle(j2,j4,j6))
			continue;

		     CHECK_EQUAL(Racah(j1,j2,j3,j4,j5,j6), Racah_NoCache(j1,j2,j3,j4,j5,j6));
		  }
	       }
	    }
	 }
      }
   }
}
