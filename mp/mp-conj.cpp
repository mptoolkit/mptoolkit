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
   if (argc != 2)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: mp-conj <psi>\n";
      return 1;
   }

   long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
   pvalue_ptr<MPWavefunction> Psi = pheap::OpenPersistent(argv[1], CacheSize);

   Psi = new MPWavefunction(conj(*Psi));

   pheap::ShutdownPersistent(Psi);
}
