// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-add2.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/matrixproduct-sum.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include <iostream>
#include "mp/copyright.h"
#include "common/environment.h"

int main(int argc, char** argv)
{
   if (argc != 5)
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "usage: mp-add2 <psi1> <psi2> <outpsi> <nstates>\n";
      return 1;
   }

   int PageSize = getenv_or_default("MP_PAGESIZE", 65536);
   long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
   pheap::Initialize(argv[3], 1, PageSize, CacheSize);
   pvalue_ptr<MPWavefunction> Psi1 = pheap::ImportHeap(argv[1]);
   pvalue_ptr<MPWavefunction> Psi2 = pheap::ImportHeap(argv[2]);
   int NStates = boost::lexical_cast<int>(argv[4]);

   std::vector<MPWavefunction> X;
   X.push_back(*Psi1);
   X.push_back(*Psi2);
   pvalue_ptr<MPWavefunction> Psi = new MPWavefunction(fancy_sum(X, NStates, NStates, 0));
   pheap::ShutdownPersistent(Psi);
}
