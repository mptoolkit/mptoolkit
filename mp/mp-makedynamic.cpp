// -*- C++ -*- $Id$

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include <iostream>
#include "mp/copyright.h"

int main(int argc, char** argv)
{
   if (argc != 6)
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "usage: mp-makedynamic <lattice> <groundstate> <frequency> <broadening> <output-operator-name>\n";
      return 1;
   }

   pvalue_ptr<OperatorList> System = pheap::OpenPersistent(argv[1], 655360);
   TRACE("here");
   pvalue_ptr<MPWavefunction> Psi = pheap::ImportHeap(argv[2]);
   double Freq = boost::lexical_cast<double>(argv[3]);
   double Broad = boost::lexical_cast<double>(argv[4]);

   MPOperator Ham = (*System)["H"];

   MPOperator Ident = (*System)["I"];

   double Energy = expectation(*Psi, Ham, *Psi).real();

   MPOperator Part = (Energy + Freq) * Ident - Ham;
   (*System.mutate())[argv[5]] = prod(Part, Part, Part.TransformsAs()) + (Broad*Broad) * Ident;

   pheap::ShutdownPersistent(System);
}
