// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-expectation-lattice.cpp
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

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpwavefunction.h"
#include "matrixproduct/mpoperatorlist.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp/copyright.h"
#include <iostream>
#include <vector>


int main(int argc, char** argv)
{
   if (argc < 5)
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "usage: mp-expectation-lattice <lattice> <psi1> <size> <operator1> [<operator2>] ... \n";
      return 1;
   }

   pvalue_ptr<OperatorList> System = pheap::OpenPersistent(argv[1], 655360, true);
   pvalue_ptr<MPWavefunction> Psi1 = pheap::ImportHeap(argv[2]);
   int L = boost::lexical_cast<int>(argv[3]);

   std::vector<std::string> Op;

   for (int i = 4; i<argc; i++)
      Op.push_back(argv[i]);


   std::cout.precision(14);

   for (int i = 1; i<=L;++i) {
      std::cout << i;

      for (unsigned j = 0; j<Op.size(); j++) {
         std::ostringstream s1;
         s1 << Op[j] << "(" << i << ')';

         double n = real(expectation(*Psi1, (*System)[s1.str()], *Psi1));

         std::cout << " " << n;
      }


      std::cout << "\n";
   }


   pheap::Shutdown();
}
