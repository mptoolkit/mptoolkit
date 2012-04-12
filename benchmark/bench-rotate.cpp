// -*- C++ -*- $Id$

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"

void DoTest(MPWavefunction& x)
{
   while (x.RightSize() > 1)
      x.RotateRight();

   while (x.LeftSize() > 1)
      x.RotateLeft();
}

int main(int argc, char** argv)
{
   if (argc != 2)
   {
      std::cerr << "usage: benchrotate <wavefunction>\n";
      return 1;
   }

   pvalue_ptr<MPWavefunction> Psi = pheap::OpenPersistent(argv[1], 655360, true);

   DoTest(*Psi.mutate());
}
