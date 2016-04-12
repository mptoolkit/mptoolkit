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
   if (argc != 3)
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "calculates the squared norm norm of |psi1> - |psi2>\n";
      std::cerr << "usage: mp-diff2 <psi1> <psi2>\n";
      return 1;
   }

   long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
   pvalue_ptr<MPWavefunction> Psi1 = pheap::OpenPersistent(argv[1], CacheSize, true);
   pvalue_ptr<MPWavefunction> Psi2 = pheap::ImportHeap(argv[2]);

   double n1 = norm_2_sq(*Psi1);
   double n2 = norm_2_sq(*Psi2);
   double ov = 2 * real(overlap(*Psi1, *Psi2));

   std::cout.precision(14);
   // if Psi1 and Psi2 are extremely close, the sum n1+n2-ov can
   // be slightly negative.
   std::cout << std::max(n1+n2-ov, 0.0) << '\n';

   pheap::Shutdown();
}
