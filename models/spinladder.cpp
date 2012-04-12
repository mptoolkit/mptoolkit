// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "models/spin.h"

std::string Coord(int i, int j)
{
   std::ostringstream S;
   S << i << ',' << j;
   return S.str();
}

int main(int argc, char** argv)
{
   if (argc != 5)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: spinladder-u1 <Spin> <L> <N-legs> <outfile>\n"
                << "Spin = spin of the bulk (half-integer)\n"
                << "L = length of ladder (number of rungs)\n"
                << "N-legs = number of legs\n"
                << "outfile = file name for output lattice.\n";
      return 1;
   }

   half_int S = boost::lexical_cast<half_int>(argv[1]);
   int L = boost::lexical_cast<int>(argv[2]);
   int N = boost::lexical_cast<int>(argv[3]);
   std::string OutFile = argv[4];

   // Construct the site block
   SiteBlock Site = CreateSpinSite(S);

   // construct the lattice
   Lattice Leg = repeat(Site, N);
   Leg.fix_coordinates();

   Lattice MyLattice = repeat(Leg, L);
   MyLattice.fix_coordinates_prepend();

   // construct the operator list for the lattice
   OperatorList OpList(MyLattice);

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number
   OperatorAtSite<OperatorList const> Sp(OpList, "Sp");
   OperatorAtSite<OperatorList const> Sm(OpList, "Sm");
   OperatorAtSite<OperatorList const> Sz(OpList, "Sz");
   OperatorAtSite<OperatorList const> Sx(OpList, "Sx");
   OperatorAtSite<OperatorList const> Sy(OpList, "Sy");
   MPOperator& H_xx_l = OpList["H_xx_l"];
   MPOperator& H_xx_t = OpList["H_xx_t"];
   MPOperator& H_yy_l = OpList["H_yy_l"];
   MPOperator& H_yy_t = OpList["H_yy_t"];
   MPOperator& H_zz_l = OpList["H_zz_l"];
   MPOperator& H_zz_t = OpList["H_zz_t"];
   MPOperator& TotalSx = OpList["Sx"];
   MPOperator& TotalSy = OpList["Sy"];
   MPOperator& TotalSz = OpList["Sz"];

   // longitudinal matrix elements
   for (int i = 1; i < L; ++i)
   {
      for (int j = 1; j <= N; ++j)
      {
	 H_xx_l += Sx(i,j) * Sx(i+1,j);
	 H_yy_l += Sy(i,j) * Sy(i+1,j);
	 H_zz_l += Sz(i,j) * Sz(i+1,j);
      }
   }

   // transverse matrix elements
   for (int i = 1; i <= L; ++i)
   {
      for (int j = 1; j < N; ++j)
      {
	 H_xx_t += Sx(i,j) * Sx(i,j+1);
	 H_yy_t += Sy(i,j) * Sy(i,j+1);
	 H_zz_t += Sz(i,j) * Sz(i,j+1);
      }
   }

   for (int i = 1; i <= L; ++i)
   {
      for (int j = 1; j <= N; ++j)
      {
	 TotalSx += Sx(i,j);
	 TotalSy += Sy(i,j);
	 TotalSz += Sz(i,j);
      }
   }

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);

   pheap::Initialize(OutFile, 1, 65536, 655360);
   pheap::ShutdownPersistent(OList);
}
