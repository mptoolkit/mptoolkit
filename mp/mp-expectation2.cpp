// -*- C++ -*- $Id$

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpwavefunction.h"
#include "matrixproduct/mpoperatorlist.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp/copyright.h"
#include <iostream>
#include "common/environment.h"

int main(int argc, char** argv)
{
   if (argc < 5 || argc > 6)
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "usage: mp-expectation2 <lattice> <psi1> <operator1> <operator2> [<psi2>]\n";
      return 1;
   }

   long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
   pvalue_ptr<OperatorList> System = pheap::OpenPersistent(argv[1], CacheSize, true);
   pvalue_ptr<MPWavefunction> Psi1 = pheap::ImportHeap(argv[2]);
   std::string Operator1 = argv[3];
   std::string Operator2 = argv[4];
   pvalue_ptr<MPWavefunction> Psi2 = argc == 6 ? pheap::ImportHeap(argv[5]) : Psi1;
   QuantumNumber Ident = QuantumNumber(System->GetSymmetryList());


   MPOperator Op1 = (*System)[Operator1];
   MPOperator Op2 = (*System)[Operator2];

   MPOperator X = prod(Op1, Op2, Ident);

   std::cout.precision(14);
   std::cout << expectation(*Psi1, X, *Psi2) << '\n';

   pheap::Shutdown();
}
