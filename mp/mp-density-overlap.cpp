// -*- C++ -*- $Id$

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpwavefunction.h"
#include "matrixproduct/mpoperatorlist.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp/copyright.h"
#include <iostream>
#include "common/environment.h"

void
dm_overlap(LinearWavefunction& Psi1, LinearWavefunction& Psi2)
{
   CHECK_EQUAL(Psi1.TransformsAs(), Psi2.TransformsAs());
   LinearWavefunction::iterator i1 = Psi1.begin();
   LinearWavefunction::iterator i2 = Psi2.begin();
   MatrixOperator M1 = MatrixOperator::make_identity(i1->Basis1());
   MatrixOperator M2 = MatrixOperator::make_identity(i2->Basis1());
   MatrixOperator E = MatrixOperator::make_identity(i2->Basis1());
   // normalize E
   E *= (1.0 / std::sqrt(trace(E)));
   std::cout << norm_frob_sq(E) << '\n';
   while (i1 != Psi1.end())
   {
      *i1 = prod(M1, *i1);
      *i2 = prod(M2, *i2);
      M1 = TruncateBasis2(*i1);
      M2 = TruncateBasis2(*i2);
      E = operator_prod(herm(*i1), E, *i2);
      double x = norm_frob_sq(triple_prod(herm(M1), E, M2));
      std::cout << x << '\n';
      ++i1;
      ++i2;
   }
}

int main(int argc, char** argv)
{
   if (argc != 3)
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "usage: mp-density-overlap <psi1> <psi2>\n";
      return 1;
   }

   long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
   pvalue_ptr<MPWavefunction> Psi1 = pheap::OpenPersistent(argv[1], CacheSize, true);
   pvalue_ptr<MPWavefunction> Psi2 = pheap::ImportHeap(argv[2]);

   std::cout.precision(14);
   MPWavefunction p1 = *Psi1;
   MPWavefunction p2 = *Psi2;
   dm_overlap(p1, p2);

   pheap::Shutdown();
}
