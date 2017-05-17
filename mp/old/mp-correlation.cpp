// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-correlation.cpp
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
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/mpwavefunction.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include <iostream>
#include <iomanip>
#include "mp/copyright.h"
#include "common/environment.h"

int main(int argc, char** argv)
{
   if (argc < 7 || argc > 8)
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "usage: mp-correlation <lattice> <psi1> <operator1> <first> <operator2> <last> [<psi2>]\n";
      return 1;
   }

   long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
   pvalue_ptr<OperatorList> System = pheap::OpenPersistent(argv[1], CacheSize, true);
   pvalue_ptr<MPWavefunction> Psi1 = pheap::ImportHeap(argv[2]);
   std::string Op1 = argv[3];
   int FirstSite = boost::lexical_cast<int>(argv[4]);
   std::string Op2 = argv[5];
   int LastSite = boost::lexical_cast<int>(argv[6]);
   pvalue_ptr<MPWavefunction> Psi2 = argc == 8 ? pheap::ImportHeap(argv[7]) : Psi1;

   // we would like to multiply the left-hand wavefunction by op1 and then calculate
   // the expectation value (Ps1Op1, Op2, Psi2), but we cannot, because
   // Op1 might not be hermitian, and we don't know how to take the adjoint of
   // an MPOperator yet...
   std::cout.precision(12);
   //std::cout.setf(std::ios_base::showpoint | std::ios_base::scientific);
   for (int Pos = FirstSite; Pos <= LastSite; ++Pos)
   {
      MPOperator Op = prod(System->Lookup(Op1, boost::lexical_cast<std::string>(FirstSite)),
                           System->Lookup(Op2, boost::lexical_cast<std::string>(Pos)),
                           QuantumNumber(System->GetSymmetryList()));
      std::complex<double> x = expectation(*Psi1, Op, *Psi2);
      std::cout << std::setw(4) << FirstSite << ' ' << std::setw(4) << Pos
                << ' ' << std::setw(18) << x.real() << ' ' << std::setw(18) << x.imag() << '\n';
   }

   pheap::Shutdown();
}
