// -*- C++ -*- $Id$

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
      print_copyright(std::cerr);
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
