// -*- C++ -*- $Id$

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/mpwavefunction.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include <iostream>
#include "mp/copyright.h"
#include "common/environment.h"

int main(int argc, char** argv)
{
   if (argc != 4)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: mp-sub <psi1> <psi2> <outpsi>\n";
      return 1;
   }

   int PageSize = getenv_or_default("MP_PAGESIZE", 65536);
   long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
   pheap::Initialize(argv[3], 1, PageSize, CacheSize);
   pvalue_ptr<MPWavefunction> Psi1 = pheap::ImportHeap(argv[1]);
   pvalue_ptr<MPWavefunction> Psi2 = pheap::ImportHeap(argv[2]);

   pvalue_ptr<MPWavefunction> Psi = new MPWavefunction(*Psi1-*Psi2);
   pheap::ShutdownPersistent(Psi);
}
