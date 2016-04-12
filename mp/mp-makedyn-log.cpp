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
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "usage: mp-makedyn-log <lattice> <groundstate-energy> <frequency> <broadening>\n";
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

   std::string OutName = std::string("GL_H(") + argv[3] + "," + argv[4] + ")";
   MPOperator Imag = Broad * Energy * Ident - Broad * Ham;
   MPOperator Real = (Energy + Freq) * Ident - Ham;
   (*System.mutate())[OutName] = Real + std::complex<double>(0.0, 1.0) * Imag;

   OutName = std::string("GL2_H(") + argv[3] + "," + argv[4] + ")";
   (*System.mutate())[OutName] = prod(Real, Real, Real.TransformsAs()) + 
      prod(Imag, Imag, Imag.TransformsAs());

   pheap::ShutdownPersistent(System);
}
