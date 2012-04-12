// -*- C++ -*- $Id$

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "pheap/pheap.h"
#include "interface/inittemp.h"

int main(int argc, char** argv)
{
   if (argc != 4)
   {
      std::cerr << "usage: make-vb-su2 <lattice> <num-sites-in-omega> <output-operator-name>\n";
      return 1;
   }

   // load the lattice file
   pvalue_ptr<OperatorList> System = pheap::OpenPersistent(argv[1], mp_pheap::CacheSize());
   int NumSitesInOmega = boost::lexical_cast<int>(argv[2]);
   std::string OpName = argv[3];

   // get the lattice size
   int LatticeSize = System->size();

   // a shortcut to refer to the spin exchange operators
   OperatorAtSite<OperatorList const, int> S(*System, "S");

   // here, Left denotes the set in Omega, Right denotes \bar\Omega 
   MPOperator SALeft, SARight, SBLeft, SBRight;
   int nAL=0, nBL=0, nAR=0, nBR=0;
   for (int i = 1; i <= NumSitesInOmega; ++i)
   {
      if (i % 2 == 0)
      {
	 ++nAL;
	 SALeft += S(i);
      }
      else
      {
	 ++nBL;
	 SBLeft += S(i);
      }
   }
   for (int i = NumSitesInOmega+1; i <= LatticeSize; ++i)
   {
      if (i % 2 == 0)
      {
	 ++nAR;
	 SARight += S(i);
      }
      else
      {
	 ++nBR;
	 SBRight += S(i);
      }
   }

   QuantumNumber Ident = QuantumNumber(System->GetSymmetryList());

   // insert our operator into the lattice
   (*System.mutate())[OpName] = sqrt(3.0) * (prod(SALeft, SBRight, Ident)
					     + prod(SBLeft, SARight, Ident));
   //      + 0.25 * (nAL*nBR + nBL*nAR) * (*System)["I"];

   // save the lattice
   pheap::ShutdownPersistent(System);
}
