// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "models/tj-u1u1.h"

typedef std::complex<double> complex;

int main(int argc, char** argv)
{
   if (argc != 3)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: tj-zigzag-u1u1 <L> <outfile>\n"
                << "L = number of lattice sites\n"
                << "outfile = file name for output lattice.\n";
      return 1;
   }

   int L = boost::lexical_cast<int>(argv[1]);
   std::string OutFile = argv[2];

   TRACE(L);

   // Construct the site block
   SiteBlock Site = CreateU1tJSite();

   // construct a lattice of L copies of Site
   Lattice MyLattice = repeat(Site, L);
   MyLattice.fix_coordinates();

   // construct the operator list for the lattice
   OperatorList OpList(MyLattice);

   OperatorAtSite<OperatorList const, int> CHup(OpList, "CHup");
   OperatorAtSite<OperatorList const, int> CHdown(OpList, "CHdown");
   OperatorAtSite<OperatorList const, int> Cup(OpList, "Cup");
   OperatorAtSite<OperatorList const, int> Cdown(OpList, "Cdown");
   OperatorAtSite<OperatorList const, int> P(OpList, "P");
   OperatorAtSite<OperatorList const, int> N(OpList, "N");
   OperatorAtSite<OperatorList const, int> Sp(OpList, "Sp");
   OperatorAtSite<OperatorList const, int> Sm(OpList, "Sm");
   OperatorAtSite<OperatorList const, int> Sz(OpList, "Sz");
   OperatorAtSite<OperatorList const, int> LocalPg(OpList, "Pg");
   MPOperator& Hamiltonian = OpList["H"];
   MPOperator& H_t_u = OpList["H_t_u"];
   MPOperator& H_t_d = OpList["H_t_d"];
   MPOperator& H_t = OpList["H_t"];
   MPOperator& H_tt_u = OpList["H_tt_u"];
   MPOperator& H_tt_d = OpList["H_tt_d"];
   MPOperator& H_tt = OpList["H_tt"];
   MPOperator& N_1 = OpList["N_1"];
   MPOperator& N_2 = OpList["N_2"];
   MPOperator& H_d = OpList["H_d"];
   MPOperator& TotalSp = OpList["Sp"];
   MPOperator& TotalSm = OpList["Sm"];
   MPOperator& TotalSz = OpList["Sz"];
   MPOperator& TotalS2 = OpList["S2"];

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number
   // hopping matrix elements

   for (int i = 1; i < L; ++i)
   {
      H_t_u -= prod(CHup(i), Cup(i+1), Ident) - prod(Cup(i), CHup(i+1), Ident);
      H_t_d -= prod(CHdown(i), Cdown(i+1), Ident) - prod(Cdown(i), CHdown(i+1), Ident);
      std::cout << "Working.... " << i << "\n";
   }

   for (int i = 1; i < L-1; ++i)
   {
      H_tt_u -= prod(CHup(i), Cup(i+2), Ident) - prod(Cup(i), CHup(i+2), Ident);
      H_tt_d -= prod(CHdown(i), Cdown(i+2), Ident) - prod(Cdown(i), CHdown(i+2), Ident);
      std::cout << "Working.... " << i << "\n";
   }

   for (int i = 1; i <= L; ++i)
   {
      if (i%2 == 1)
	 N_1 += N(i);
      else
	 N_2 += N(i);

      TotalSp += Sp(i);
      TotalSm += Sm(i);
      TotalSz += Sz(i);
   }

   H_t = H_t_u + H_t_d;
   H_tt = H_tt_u + H_tt_d;

   H_d = 0.5 * (N_1 - N_2);

   TotalS2 = 0.5*(TotalSp*TotalSm + TotalSm*TotalSp) + TotalSz*TotalSz;

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);

   pheap::Initialize(OutFile, 1, 65536, 655360);
   pheap::ShutdownPersistent(OList);
}
