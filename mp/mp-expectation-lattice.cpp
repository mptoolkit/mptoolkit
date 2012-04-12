// -*- C++ -*- $Id$

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpwavefunction.h"
#include "matrixproduct/mpoperatorlist.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp/copyright.h"
#include <iostream>
#include <vector>


int main(int argc, char** argv)
{
   if (argc < 5)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: mp-expectation-lattice <lattice> <psi1> <size> <operator1> [<operator2>] ... \n";
      return 1;
   }

   pvalue_ptr<OperatorList> System = pheap::OpenPersistent(argv[1], 655360, true);
   pvalue_ptr<MPWavefunction> Psi1 = pheap::ImportHeap(argv[2]);
   int L = boost::lexical_cast<int>(argv[3]);
   
   std::vector<std::string> Op;
   
   for (int i = 4; i<argc; i++)
      Op.push_back(argv[i]);

   
   std::cout.precision(14);
   
   for (int i = 1; i<=L;++i) {
      std::cout << i;
      
      for (unsigned j = 0; j<Op.size(); j++) {
         std::ostringstream s1;
         s1 << Op[j] << "(" << i << ')';

         double n = real(expectation(*Psi1, (*System)[s1.str()], *Psi1));
         
         std::cout << " " << n;
      }
         
      
      std::cout << "\n";
   }
   

   pheap::Shutdown();
}
