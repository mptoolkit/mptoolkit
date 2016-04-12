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
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "usage: mp-defrag <psi>\n";
      return 1;
   }

   long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
   pvalue_ptr<MPWavefunction> Psi = pheap::OpenPersistent(argv[1], CacheSize);
   {
      pvalue_ptr<MPWavefunction>::lock_type l(Psi);
      // Make a copy of the handles, to keep the originals on disk while we make a copy
      std::vector<pvalue_handle<MPWavefunction> > Keep(l->base_begin(), l->base_end());
      for (MPWavefunction::base_iterator I = l->base_begin(); I != l->base_end(); ++I)
         I->lock().mutate();
   }
   pheap::ShutdownPersistent(Psi);
}
