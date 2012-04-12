// -*- C++ -*- $Id$

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include <iostream>
#include "mp/copyright.h"

int main(int argc, char** argv)
{
   if (argc != 5)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: mp-makedyn <lattice> <groundstate-energy> <frequency> <broadening>\n";
      return 1;
   }

   pvalue_ptr<OperatorList> System = pheap::OpenPersistent(argv[1], 655360);
   TRACE("here");
   //   pvalue_ptr<MPWavefunction> Psi = pheap::ImportHeap(argv[2]);
   double Freq = boost::lexical_cast<double>(argv[3]);
   double Broad = boost::lexical_cast<double>(argv[4]);

   MPOperator Ham = (*System)["H"];

   MPOperator Ident = (*System)["I"];

   double Energy = boost::lexical_cast<double>(argv[2]); // expectation(*Psi, Ham, *Psi).real();

   std::string OutName = std::string("G_H(") + argv[3] + "," + argv[4] + ")";
   (*System.mutate())[OutName] = std::complex<double>(Energy + Freq, Broad) * Ident - Ham;

   OutName = std::string("G2_H(") + argv[3] + "," + argv[4] + ")";
   MPOperator Part = (Energy + Freq) * Ident - Ham;
   (*System.mutate())[OutName] = prod(Part, Part, Part.TransformsAs()) + (Broad*Broad) * Ident;

   pheap::ShutdownPersistent(System);
}
