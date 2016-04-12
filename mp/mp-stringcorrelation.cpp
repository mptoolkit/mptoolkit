// -*- C++ -*- $Id$

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/mpwavefunction.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include <iostream>
#include <iomanip>
#include "mp/copyright.h"
#include "common/environment.h"
#include "interface/inittemp.h"

int main(int argc, char** argv)
{
   if (argc < 8 || argc > 10)
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "usage: mp-stringcorrelation <lattice> <psi1> <operator1> <first> <join> <operator2> <last> [<psi2>] [joinfirst]\n";
      std::cerr << "if joinfirst exists, then the join operator will be applied to the first site.\n";
      return 1;
   }


   mp_pheap::InitializeTempPHeap(false);
   pvalue_ptr<OperatorList> System = pheap::ImportHeap(argv[1]);
   pvalue_ptr<MPWavefunction> Psi1 = pheap::ImportHeap(argv[2]);
   std::string Op1 = argv[3];
   int FirstSite = boost::lexical_cast<int>(argv[4]);
   std::string Join = argv[5];
   std::string Op2 = argv[6];
   int LastSite = boost::lexical_cast<int>(argv[7]);
   pvalue_ptr<MPWavefunction> Psi2 = argc >= 9 ? pheap::ImportHeap(argv[8]) : Psi1;
   bool JoinFirst = 0;
   if (argc == 10)
      JoinFirst = true;

   // we would like to multiply the left-hand wavefunction by op1 and then calculate
   // the expectation value (Ps1Op1, Op2, Psi2), but we cannot, because 
   // Op1 might not be hermitian, and we don't know how to take the adjoint of
   // an MPOperator yet...
   std::cout.precision(14);
   MPOperator Operator1 = (Op1 == "I") ? (*System)["I"] 
      : System->Lookup(Op1, boost::lexical_cast<std::string>(FirstSite));
   if (JoinFirst)
   {
      // incorporate the join into the first site as well
      Operator1 = prod(Operator1, System->Lookup(Join, boost::lexical_cast<std::string>(FirstSite)), 
                       Operator1.TransformsAs());
   }

   for (int Pos = FirstSite+1; Pos <= LastSite; ++Pos)
   {
      MPOperator Operator2 = (Op2 == "I") ? (*System)["I"]
         : System->Lookup(Op2, boost::lexical_cast<std::string>(Pos));
      MPOperator Op = prod(Operator1, Operator2, QuantumNumber(System->GetSymmetryList()));
      for (int P = FirstSite+1; P < Pos; ++P)
      {
	 Op = prod(Op, System->Lookup(Join, boost::lexical_cast<std::string>(P)), 
		   QuantumNumber(System->GetSymmetryList()));
      }
      std::cout << std::setw(4) << FirstSite << ' ' << std::setw(4) << Pos 
		<< ' ' << std::setw(16)
		<< expectation(*Psi1, Op, *Psi2) << '\n';
   }

   pheap::Shutdown();
}
