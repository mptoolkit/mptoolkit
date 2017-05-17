// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-rme.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
#include "matrixproduct/mpwavefunction.h"
#include "matrixproduct/mpoperatorlist.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp/copyright.h"
#include <iostream>
#include "common/environment.h"
#include "interface/operator-parser.h"

int main(int argc, char** argv)
{
   if (argc < 3 || argc > 4)
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "usage: mp-rme <psi1> <operator> [<psi2>]\n";
      return 1;
   }

   long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
   pvalue_ptr<MPWavefunction> Psi1 = pheap::OpenPersistent(argv[1], CacheSize, true);
   MPOperator Op = ParseOperator(argv[2]);
   pvalue_ptr<MPWavefunction> Psi2 = argc == 4 ? pheap::ImportHeap(argv[3]) : Psi1;

   std::cout.precision(14);
   MPStateComponent x = reduced_matrix_element(*Psi1, Op, *Psi2);
   std::cout << trace(x[0](0,0)) << '\n';

   pheap::Shutdown();
}
