// -*- C++ -*- $Id$

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "quantumnumbers/u1.h"
#include "pheap/pheap.h"
#include "matrixproduct/copyright.h"
#include <iostream>

int main(int argc, char** argv)
{
   if (argc != 6)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: mp-itime <lattice> <operator> <psi> <delta> <iterations>\n";
      return 1;
   }

   pvalue_ptr<MPWavefunction> Psi = pheap::OpenPersistent(argv[3], 655360);
   pvalue_ptr<OperatorList> System = pheap::ImportHeap(argv[1]);
   std::string Operator = argv[2];
   double Delta = boost::lexical_cast<double>(argv[4]);
   double Iterations = boost::lexical_cast<int>(argv[5]);

   MPWavefunction P = *Psi;
   P.normalize();
   MPOperator H = (*System)[Operator];

   MPWavefunction HP = prod(H, P, P.TransformsAs());

   std::cout.precision(14);
   std::cout << "Iteration                 Energy                Overlap\n";
   for (int i = 0; i < Iterations; ++i)
   {
      std::cout << std::setw(9) << i+1 << "   "
		<< std::setw(20) << overlap(P, HP);
      HP *= Delta;
      P = P + HP;
      P.normalize();
      std::cout << "   " << std::setw(20) << overlap(*Psi, P) << '\n';
      HP = prod(H, P, P.TransformsAs());
   }

   *Psi.mutate() = P;
   pheap::ShutdownPersistent(Psi);
}
