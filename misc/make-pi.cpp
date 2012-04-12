// -*- C++ -*- $Id$

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "pheap/pheap.h"

int main(int argc, char** argv)
{
   if (argc != 2)
   {
      std::cerr << "usage: make-pi <lattice>\n";
      return 1;
   }

   // load the lattice file
   pvalue_ptr<OperatorList> System = pheap::OpenPersistent(argv[1], 655360);

   // get the lattice size
   int LatticeSize = System->size();

   // a shortcut to refer to the "S" (spin) operator   
   OperatorAtSite<OperatorList const, int> S(*System, "S");

   // construct our output operator
   MPOperator Spi;

   for (int i = 1; i <= LatticeSize; ++i)
   {
      Spi += pow(-1, i) * S(i);
   }

   // insert our operator into the lattice
   (*System.mutate())["Spi"] = Spi;

   // save the lattice
   pheap::ShutdownPersistent(System);
}
