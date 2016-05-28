// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-trotter.cpp
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
#include "matrixproduct/mpexponential.h"
#include "matrixproduct/operatoratsite.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include <iostream>
#include <iomanip>
#include "mp/copyright.h"

int main(int argc, char** argv)
{
   if (argc < 3 || argc > 3)
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "usage: mp-trotter <lattice> <x>\n"
                << "Calculates TrotterE(x) and TrotterO(x) for even and odd Trotter slices.\n";
      std::cerr << "Requires the lattice contains operators Bond(i)\n";
      return 1;
   }

   pvalue_ptr<OperatorList> System = pheap::OpenPersistent(argv[1], 655360);
   std::complex<double> x = boost::lexical_cast<std::complex<double> >(argv[2]);

   std::string tEStr = std::string("TrotterE(") + argv[2] + ")";
   std::string tOStr = std::string("TrotterO(") + argv[2] + ")";

   OperatorAtSite<OperatorList const, int> B(*System, "Bond");
   MPOperator tE = (*System)["I"];
   MPOperator tO = (*System)["I"];

   QuantumNumber Ident = tE.TransformsAs();
   int i = 0;
   while (i < System->size() && !B.HasOperator(i)) ++i;
   CHECK(B.HasOperator(i))(i);
   bool Found = true;
   while (Found)
   {
      if (B.HasOperator(i))
      {
         if (i%2 == 0)
            tE = prod(tE, BondExponential(x, SplitOperator(B(i))).AsMPOperator(), Ident);
         else
            tO = prod(tO, BondExponential(x, SplitOperator(B(i))).AsMPOperator(), Ident);
      }
      else Found = false;
      ++i;
      //      Found = false;
   }

   (*System.mutate())[tEStr] = tE;
   (*System.mutate())[tOStr] = tO;

   std::cout << "Created operator: " << tEStr << "\nCreated operator: " << tOStr << std::endl;

   pheap::ShutdownPersistent(System);
} 
