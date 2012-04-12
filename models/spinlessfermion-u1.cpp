// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "models/spinlessfermion-u1.h"

int main(int argc, char** argv)
{
   if (argc != 3)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: spinlessfermion-u1 <L> <outfile>\n"
                << "L = number of lattice sites\n"
                << "outfile = file name for output lattice.\n";
      std::cerr << "\nOperators:\n"
                << "H_t = nearest-neighbor hopping\n"
                << "H_tt = next-nearest-neighbor hopping\n"
                << "H_V = nearest-neighbor coulomb repulsion\n"
                << "N_1 = particle number of odd sites N(1)+N(3)+N(5)+...\n"
                << "N_2 = particle number of even sites N(2)+N(4)+N(6)+...\n"
                << "N_d = 0.5*(N_1 - N_2)\n"
         ;
      return 1;
   }

   int L = boost::lexical_cast<int>(argv[1]);
   std::string OutFile(argv[2]);
   
   // Construct the site block
   SiteBlock Site = CreateU1SpinlessFermion();

   // construct a lattice of L copies of Site
   Lattice MyLattice = repeat(Site, L);
   MyLattice.fix_coordinates();

   // construct the operator list for the lattice
   OperatorList OpList(MyLattice);
   OperatorAtSite<OperatorList const> CH(OpList, "CH");
   OperatorAtSite<OperatorList const> C(OpList, "C");
   OperatorAtSite<OperatorList const> N(OpList, "N");
   MPOperator& H_t = OpList["H_t"];
   MPOperator& H_tt = OpList["H_tt"];
   MPOperator& H_V = OpList["H_V"];
   MPOperator& N_1 = OpList["N_1"];
   MPOperator& N_2 = OpList["N_2"];
   MPOperator& H_d = OpList["H_d"];

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number
   // interaction matrix elements
   for (int i = 1; i <= L-1; ++i)
   {
      H_t -= prod(CH(i), C(i+1), Ident) + prod(CH(i+1), C(i), Ident);
      H_V += N(i) * N(i+1);
      std::cout << "Working.... " << i << "\n";
   }
   for (int i = 1; i <= L-2; ++i)
   {
      H_tt -= prod(CH(i), C(i+2), Ident) + prod(CH(i+2), C(i), Ident);
      std::cout << "Working.... " << i << "\n";
   }
   for (int i = 1; i <= L; ++i)
   {
      if (i%2 == 1)
	 N_1 += N(i);
      else
	 N_2 += N(i);
   }

   H_d = 0.5 * (N_1 - N_2);

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);

   pheap::Initialize(OutFile, 1, 65536, 655360);
   pheap::ShutdownPersistent(OList);
}
