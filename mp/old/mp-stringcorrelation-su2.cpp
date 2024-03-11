// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/old/mp-stringcorrelation-su2.cpp
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
#include "quantumnumbers/all_symmetries.h"
#include "matrixproduct/operatoratsite.h"
#include "pheap/pheap.h"
#include <iostream>
#include <iomanip>
#include "mp/copyright.h"

// project out the S=2 bond
MPOperator Projector(OperatorList const& System, int i1, int i2)
{
   OperatorAtSite<OperatorList const, int> S(System, "S");
   MPOperator Bond = -sqrt(3.0) * prod(S(i1), S(i2), QuantumNumber(System.GetSymmetryList()));
   return (1.0/3.0) * System["I"] - Bond - (1.0/3.0)*prod(Bond,Bond,Bond.TransformsAs());
}

// project onto the singlet bond
MPOperator ProjectorSinglet(OperatorList const& System, int i1, int i2)
{
   OperatorAtSite<OperatorList const, int> S(System, "S");
   MPOperator Bond = -sqrt(3.0) * prod(S(i1), S(i2), QuantumNumber(System.GetSymmetryList()));
   return (-1.0/3.0) * (System["I"] + Bond);
}

// project onto the triplet bond
MPOperator ProjectorTriplet(OperatorList const& System, int i1, int i2)
{
   OperatorAtSite<OperatorList const, int> S(System, "S");
   MPOperator Bond = -sqrt(3.0) * prod(S(i1), S(i2), QuantumNumber(System.GetSymmetryList()));
   return System["I"] - 0.5 * (Bond + prod(Bond, Bond, Bond.TransformsAs()));
}

// project onto the quintuplet bond
MPOperator ProjectorQuintuplet(OperatorList const& System, int i1, int i2)
{
   OperatorAtSite<OperatorList const, int> S(System, "S");
   MPOperator Bond = -sqrt(3.0) * prod(S(i1), S(i2), QuantumNumber(System.GetSymmetryList()));
   return (1.0/3.0)*System["I"] + 0.5*Bond + (1.0/6.0)*prod(Bond, Bond, Bond.TransformsAs());
}

MPOperator BondPhase(OperatorList const& System, int i1, int i2)
{
   OperatorAtSite<OperatorList const, int> S(System, "S");
   MPOperator Bond = -sqrt(3.0) * prod(S(i1), S(i2), QuantumNumber(System.GetSymmetryList()));
   return -1.0*System["I"] + Bond + prod(Bond,Bond,Bond.TransformsAs());
}

int main(int argc, char** argv)
{
   if (argc < 5 || argc > 6)
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "usage: mp-stringcorrelation-su2 <lattice> <psi1> <first> <last> [<psi2>]\n";
      return 1;
   }

   pvalue_ptr<OperatorList> System = pheap::OpenPersistent(argv[1], 655360, true);
   pvalue_ptr<MPWavefunction> Psi1 = pheap::ImportHeap(argv[2]);
   int FirstSite = boost::lexical_cast<int>(argv[3]);
   int LastSite = boost::lexical_cast<int>(argv[4]);
   pvalue_ptr<MPWavefunction> Psi2 = argc == 6 ? pheap::ImportHeap(argv[5]) : Psi1;

#if 0
   MPOperator P = (*System)["I"];
   for (int Pos = FirstSite; Pos < LastSite; ++Pos)
   {
      P = prod(P, Projector(*System, Pos, Pos+1), P.TransformsAs());
      std::cout << std::setw(4) << FirstSite << ' ' << std::setw(4) << (Pos+1)
                << ' ' << std::setw(16)
                << expectation(*Psi1, P, *Psi2) << '\n';
   }
#endif

#if 0
   OperatorAtSite<OperatorList const, int> S(*System, "S");
   MPOperator I = (*System)["I"];
   QuantumNumber Ident = QuantumNumber(System->GetSymmetryList());
   std::cout.precision(14);
   MPOperator Op = S(FirstSite);
   for (int Pos = FirstSite+1; Pos <= LastSite; ++Pos)
   {
      MPOperator ThisOp = prod(Op, S(Pos), Ident);
      std::cout << std::setw(4) << FirstSite << ' ' << std::setw(4) << Pos
                << ' ' << std::setw(16)
                << (-std::sqrt(3.0) * expectation(*Psi1, ThisOp, *Psi2)) << std::endl;

      Op = prod(Op, BondPhase(*System, Pos-1, Pos), Op.TransformsAs());
   }
#endif

#if 1
   OperatorAtSite<OperatorList const, int> S(*System, "S");
   MPOperator I = (*System)["I"];
   QuantumNumber Ident = QuantumNumber(System->GetSymmetryList());

   std::complex<double> x(0.0, 0.0);
   std::cout.precision(14);
   MPOperator Op = S(FirstSite);
   for (int Pos = FirstSite+1; Pos <= LastSite; ++Pos)
   {
      MPOperator ThisOp = prod(Op, S(Pos), Ident);
      std::cout << std::setw(4) << FirstSite << ' ' << std::setw(4) << Pos
                << ' ' << x << ' ' << std::setw(16)
                << (-std::sqrt(3.0) * expectation(*Psi1, ThisOp, *Psi2)) << '\n';

      //      Op += x / fabs(x) * std::sqrt(fabs(x))
      //         * prod(Op, prod(S(Pos), S(Pos), QuantumNumber(Op.GetSymmetryList(), "2")),
      //                Op.TransformsAs());

      Op += x * prod(Op, S(Pos), Op.TransformsAs());

      //Op = prod(Op, prod(S(Pos), S(Pos), QuantumNumber(Op.GetSymmetryList(), "2")),
      //          Op.TransformsAs());
      //Op *= std::sqrt(27.0 / 20.0);
      //Op = prod(Op, I, Op.TransformsAs());
   }
#endif

   pheap::Shutdown();
}
