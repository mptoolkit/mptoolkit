// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "models/spin.h"

int main(int argc, char** argv)
{
   if (argc != 5)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: spinchain <L> <EdgeSpin> <BulkSpin> <outfile>\n"
                << "L = number of lattice sites\n"
                << "EdgeSpin = spin of the left and right extreme sites\n"
                << "BulkSpin = spin away from the boundary\n"
                << "outfile = file name for output lattice.\n";
      return 1;
   }

   int L = boost::lexical_cast<int>(argv[1]);
   half_int Edge = boost::lexical_cast<half_int>(argv[2]);
   half_int Bulk = boost::lexical_cast<half_int>(argv[3]);
   std::string OutName = argv[4];

   // Construct the site block
   SiteBlock EdgeSite = CreateSpinSite(Edge);
   SiteBlock BulkSite = CreateSpinSite(Bulk);

   // construct a lattice of L copies of Site
   Lattice MyLattice = join(EdgeSite, repeat(BulkSite, L-2), EdgeSite);
   MyLattice.fix_coordinates();

   // construct the operator list for the lattice
   OperatorList OpList(MyLattice);

   OperatorAtSite<OperatorList const, int> Sp(OpList, "Sp");
   OperatorAtSite<OperatorList const, int> Sm(OpList, "Sm");
   OperatorAtSite<OperatorList const, int> Sz(OpList, "Sz");
   OperatorAtSite<OperatorList const, int> Sx(OpList, "Sx");
   OperatorAtSite<OperatorList const, int> Sy(OpList, "Sy");
   MPOperator& H_xx = OpList["H_xx"];
   MPOperator& H_yy = OpList["H_yy"];
   MPOperator& H_zz = OpList["H_zz"];
   MPOperator& TotalSp = OpList["Sp"];
   MPOperator& TotalSm = OpList["Sm"];
   MPOperator& TotalSz = OpList["Sz"];
   MPOperator& TotalS2 = OpList["S2"];
   MPOperator& TotalSx = OpList["Sx"];
   MPOperator& TotalSy = OpList["Sy"];

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number
   // interaction matrix elements

   // This approach seems to be buggy, so for now we revert to the old style
#if 0
   H_xx = CreateRepeatedOperator(MyLattice, "Sx", "Sx");
   H_yy = CreateRepeatedOperator(MyLattice, "Sy", "Sy");
   H_zz = CreateRepeatedOperator(MyLattice, "Sz", "Sz");

   TotalSx = CreateRepeatedOperator(MyLattice, "Sx");
   TotalSy = CreateRepeatedOperator(MyLattice, "Sy");
   TotalSz = CreateRepeatedOperator(MyLattice, "Sz");
   TotalSp = CreateRepeatedOperator(MyLattice, "Sp");
   TotalSm = CreateRepeatedOperator(MyLattice, "Sm");

#else

   for (int i = 1; i < L; ++i)
   {
      H_xx += Sx(i) * Sx(i+1);
      H_yy += Sy(i) * Sy(i+1);
      H_zz += Sz(i) * Sz(i+1);
   }

   for (int i = 1; i <= L; ++i)
   {
      TotalSx += Sx(i);
      TotalSy += Sy(i);
      TotalSz += Sz(i);
      TotalSp += Sp(i);
      TotalSm += Sm(i);
   }
#endif

   TotalS2 = 0.5 * (TotalSp * TotalSm + TotalSm * TotalSp)
      + TotalSz * TotalSz;

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);

   pheap::Initialize(OutName, 1, 65536, 655360);
   pheap::ShutdownPersistent(OList);
}
