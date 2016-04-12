// -*- C++ -*- $Id$

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include <iostream>
#include "mp/copyright.h"

int main(int argc, char** argv)
{
   if (argc < 4 || argc > 5)
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "usage: mp-expectation <lattice> <psi1> <operator> [<psi2>]\n";
      return 1;
   }

   pvalue_ptr<OperatorList> System = pheap::OpenPersistent(argv[1], 655360, true);
   pvalue_ptr<MPWavefunction> Psi1 = pheap::ImportHeap(argv[2]);
   std::string Operator = argv[3];
   pvalue_ptr<MPWavefunction> Psi2 = argc == 5 ? pheap::ImportHeap(argv[4]) : Psi1;

   std::cout.precision(14);
   std::cout << expectation_r(*Psi1, (*System)[Operator], *Psi2) << '\n';

   pheap::Shutdown();
}
