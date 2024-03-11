// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/old/mp-theta2.cpp
//
// Copyright (C) 2004-2017 Ian McCulloch <ian@qusim.net>
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

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/mpwavefunction.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include <iostream>
#include "mp/copyright.h"
#include "common/environment.h"

int main(int argc, char** argv)
{
   if (argc != 3)
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "calculates the angle theta^2 given by |<psi1|psi2>| = ||psi1|| ||psi2|| cos theta\n";
      std::cerr << "usage: mp-theta2 <psi1> <psi2>\n";
      return 1;
   }

   long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
   pvalue_ptr<MPWavefunction> Psi1 = pheap::OpenPersistent(argv[1], CacheSize, true);
   pvalue_ptr<MPWavefunction> Psi2 = pheap::ImportHeap(argv[2]);

   double n1 = norm_2_sq(*Psi1);
   double n2 = norm_2_sq(*Psi2);
   double ov = 2 * real(overlap(*Psi1, *Psi2));

   std::cout.precision(14);
   // if Psi1 and Psi2 are extremely close, the sum n1+n2-ov can
   // be slightly negative.
   std::cout << (2.0 - ov / std::sqrt(n1 * n2)) << '\n';

   pheap::Shutdown();
}
