// -*- C++ -*- $Id: hubbard-u1su2.cpp 819 2008-01-14 05:03:01Z ianmcc $

#include "pheap/pheap.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "models/tj-u1su2.h"

typedef std::complex<double> complex;

int main(int argc, char** argv)
{
   if (argc != 3)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: tj-zigzag-u1su2 <L> <outfile>\n"
                << "L = number of lattice sites\n"
                << "outfile = file name for output lattice.\n";
      return 1;
   }

   int L = boost::lexical_cast<int>(argv[1]);
   std::string OutFile = argv[2];

   TRACE(L);

   // Construct the site block
   SiteBlock Site = CreateU1SU2tJSite();

   // construct a lattice of L copies of Site
   Lattice MyLattice = repeat(Site, L);
   MyLattice.fix_coordinates();

   // construct the operator list for the lattice
   OperatorList OpList(MyLattice);

   OperatorAtSite<OperatorList const, int> CH(OpList, "CH");
   OperatorAtSite<OperatorList const, int> C(OpList, "C");
   OperatorAtSite<OperatorList const, int> P(OpList, "P");
   OperatorAtSite<OperatorList const, int> N(OpList, "N");
   OperatorAtSite<OperatorList const, int> LocalPg(OpList, "Pg");
   OperatorAtSite<OperatorList, int> Bond(OpList, "Bond");
   MPOperator& Hamiltonian = OpList["H"];
   MPOperator& H_t = OpList["H_t"];
   MPOperator& H_tt = OpList["H_tt"];
   MPOperator& N_1 = OpList["N_1"];
   MPOperator& N_2 = OpList["N_2"];
   MPOperator& H_d = OpList["H_d"];

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number
   // hopping matrix elements

   double Sqrt2 = -std::sqrt(2.0);

   for (int i = 1; i < L; ++i)
   {
      H_t += Sqrt2 * (prod(CH(i), C(i+1), Ident) + prod(C(i), CH(i+1), Ident));
      std::cout << "Working.... " << i << "\n";
   }

   for (int i = 1; i < L-1; ++i)
   {
      H_tt += Sqrt2 * (prod(CH(i), C(i+2), Ident) + prod(C(i), CH(i+2), Ident));
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
