// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/old/mp-generalizedstring.cpp
//
// Copyright (C) 2004-2017 Ian McCulloch <ian@qusim.net>
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

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include <iostream>
#include <iomanip>
#include "mp/copyright.h"
#include "common/environment.h"

typedef std::complex<double> complex;

int main(int argc, char** argv)
{
   if (argc < 9 || argc > 10)
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "usage: mp-stringcorrelation <lattice> <psi1> <operator1> <first> <kS> <kH> <operator2> <last> [<psi2>]\n";
      return 1;
   }

   long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
   pvalue_ptr<OperatorList> System = pheap::OpenPersistent(argv[1], CacheSize, true);
   pvalue_ptr<MPWavefunction> Psi1 = pheap::ImportHeap(argv[2]);
   std::string Op1 = argv[3];
   int FirstSite = boost::lexical_cast<int>(argv[4]);
   int kS = boost::lexical_cast<int>(argv[5]);
   int kH = boost::lexical_cast<int>(argv[6]);
   std::string Op2 = argv[7];
   int LastSite = boost::lexical_cast<int>(argv[8]);
   pvalue_ptr<MPWavefunction> Psi2 = argc == 10 ? pheap::ImportHeap(argv[9]) : Psi1;

   int const L = System->size();

   // we would like to multiply the left-hand wavefunction by op1 and then calculate
   // the expectation value (Ps1Op1, Op2, Psi2), but we cannot, because
   // Op1 might not be hermitian, and we don't know how to take the adjoint of
   // an MPOperator yet...
   std::cout.precision(14);
   MPOperator Operator1 = (Op1 == "I") ? (*System)["I"]
      : System->Lookup(Op1, boost::lexical_cast<std::string>(FirstSite));

   OperatorAtSite<OperatorList const> SpinonProj(*System, "N_S");
   OperatorAtSite<OperatorList const> HolonProj(*System, "N_H");

   for (int Pos = FirstSite+1; Pos <= LastSite; ++Pos)
   {
      MPOperator Operator2 = (Op2 == "I") ? (*System)["I"]
         : System->Lookup(Op2, boost::lexical_cast<std::string>(Pos));
      MPOperator Op = prod(Operator1, Operator2, QuantumNumber(System->GetSymmetryList()));
      for (int P = FirstSite+1; P < Pos; ++P)
      {
         MPOperator Join = std::exp(complex(0,math_const::pi * kS / L)) * SpinonProj(P)
            + std::exp(complex(0,math_const::pi * kH / L)) * HolonProj(P);

         Op = prod(Op, Join, QuantumNumber(System->GetSymmetryList()));
      }
      std::cout << std::setw(4) << FirstSite << ' ' << std::setw(4) << Pos
                << ' ' << std::setw(16)
                << expectation(*Psi1, Op, *Psi2) << '\n';
   }

   pheap::Shutdown();
}
