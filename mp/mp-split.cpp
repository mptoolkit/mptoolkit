// -*- C++ -*- $Id$

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/mpwavefunction.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include <iostream>
#include "mp/copyright.h"
#include "common/environment.h"

typedef std::complex<double> complex;

int main(int argc, char** argv)
{
   if (argc != 4)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: mp-split <psi> <real-part> <imag-part>\n";
      return 1;
   }

   long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
   pvalue_ptr<MPWavefunction> PsiPtr = pheap::OpenPersistent(argv[1], CacheSize, true);

   MPWavefunction Psi = *PsiPtr;
   MPWavefunction PsiBar = conj(Psi);

   {
      pvalue_ptr<MPWavefunction> PsiReal = new MPWavefunction(0.5 * (Psi + PsiBar));
      //   }
      //   {
      pvalue_ptr<MPWavefunction> PsiImag = 
         new MPWavefunction(complex(0.0,-0.5) * (Psi - PsiBar));
      pheap::ExportHeap(argv[3], PsiImag);
      pheap::ExportHeap(argv[2], PsiReal);
   }
}
