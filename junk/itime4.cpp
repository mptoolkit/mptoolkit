// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// junk/itime4.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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
#include "quantumnumbers/u1.h"
#include "pheap/pheap.h"
#include <iostream>
#include "matrixproduct/copyright.h"

int main(int argc, char** argv)
{
   if (argc != 6)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: mp-itime4 <lattice> <operator> <psi> <delta> <iterations>\n";
      return 1;
   }

   pvalue_ptr<MPWavefunction> Psi = pheap::OpenPersistent(argv[3], 655360);
   pvalue_ptr<OperatorList> System = pheap::ImportHeap(argv[1]);
   std::string Operator = argv[2];
   double Delta = boost::lexical_cast<double>(argv[4]);
   double Iterations = boost::lexical_cast<int>(argv[5]);

   MPWavefunction P = *Psi;
   P.normalize();
   MPOperator H = (*System)[Operator];

   std::cout << "Calculating 4th order taylor series expansion of exp[Delta * H]..." << std::endl;
   // calculate H^2
   MPOperator H2 = prod(H, H, H.TransformsAs());
   // calculate H^3
   MPOperator H3 = prod(H, H2, H.TransformsAs());
   // calculate H^4
   MPOperator H4 = prod(H, H3, H.TransformsAs());

   MPOperator Taylor = (*System)["I"]
      + Delta * H
      + std::pow(Delta, 2) * H2
      + std::pow(Delta, 3) * H3
      + std::pow(Delta, 4) * H4;

   std::cout << "done!  Starting calculation...\n";

   std::cout.precision(14);
   std::cout << "Iteration                 Energy                Overlap\n";
   for (int i = 0; i < Iterations; ++i)
   {
      P = prod(Taylor, P, P.TransformsAs());
      P.normalize();
      std::cout << std::setw(9) << i+1 << "   "
                << std::setw(20) << expectation(P, H, P)
                << "   " << std::setw(20) << overlap(*Psi, P) << '\n';
   }

   *Psi.mutate() = P;
   pheap::ShutdownPersistent(Psi);
}
