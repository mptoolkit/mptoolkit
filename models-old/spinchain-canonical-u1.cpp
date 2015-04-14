// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "models/spin-u1.h"

int main(int argc, char** argv)
{
   if (argc != 8)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: spinchain-canonical-u1 <L> <EdgeSpin> <BulkSpin> <Jperp> <Jz> <Beta> <outfile>\n"
                << "L = number of lattice sites\n"
                << "EdgeSpin = spin of the left and right extreme sites\n"
                << "BulkSpin = spin away from the boundary\n"
                << "Jperp = perpendicular coupling\n"
                << "Jz = Sz-Sz coupling\n"
                << "Beta = biquadratic coupling\n"
                << "outfile = file name for output lattice.\n";
      return 1;
   }

   int L = boost::lexical_cast<int>(argv[1]);
   half_int Edge = boost::lexical_cast<half_int>(argv[2]);
   half_int Bulk = boost::lexical_cast<half_int>(argv[3]);
   double Jperp = boost::lexical_cast<double>(argv[4]);
   double Jz = boost::lexical_cast<double>(argv[5]);
   double Beta = boost::lexical_cast<double>(argv[6]);

   TRACE(L)(Edge)(Bulk)(Jperp)(Jz)(Beta);

   CHECK(L >= 2);
   
   // Construct the site block
   SiteBlock EdgeSite = CreateU1SpinSite(Edge);
   SiteBlock BulkSite = CreateU1SpinSite(Bulk);
   SiteBlock AuxEdgeSite = CreateU1SpinSite(Edge, "Szaux");
   SiteBlock AuxBulkSite = CreateU1SpinSite(Bulk, "Szaux");

   Lattice EdgeUnitCell(SymmetryList("Sz:U(1),Szaux:U(1)"), EdgeSite, AuxEdgeSite);
   EdgeUnitCell.fix_coordinates("", "aux");
   Lattice BulkUnitCell(SymmetryList("Sz:U(1),Szaux:U(1)"), BulkSite, AuxBulkSite);
   BulkUnitCell.fix_coordinates("", "aux");
   
   // construct a lattice of L copies of Site
   Lattice MyLattice = join(EdgeUnitCell, repeat(BulkUnitCell, L-2), EdgeUnitCell);
   MyLattice.fix_coordinates_prepend();

   // construct the operator list for the lattice
   OperatorList OpList(MyLattice);

   OperatorAtSite<OperatorList const> Sp(OpList, "Sp");
   OperatorAtSite<OperatorList const> Sm(OpList, "Sm");
   OperatorAtSite<OperatorList const> Sz(OpList, "Sz");
   MPOperator& Hamiltonian = OpList["H"];
   MPOperator& AuxHamiltonian = OpList["AuxH"];
   MPOperator& TotalSp = OpList["Sp"];
   MPOperator& TotalSm = OpList["Sm"];
   MPOperator& TotalSz = OpList["Sz"];
   MPOperator& TotalS2 = OpList["S2"];

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number
   // interaction matrix elements
   for (int i = 1; i < L; ++i)
   {
      MPOperator SS = 0.5 * Jperp * (prod(Sp(i), Sm(i%L+1), Ident)
                                     + prod(Sm(i), Sp(i%L+1), Ident))
         + Jz * prod(Sz(i), Sz(i%L+1), Ident);
      if (Beta != 0)
      {
         MPOperator SSIso = 0.5 * (prod(Sp(i), Sm(i%L+1), Ident)
                                   + prod(Sm(i), Sp(i%L+1), Ident))
            + prod(Sz(i), Sz(i%L+1), Ident);
         SS += Beta * prod(SSIso, SSIso, Ident);
      }

      Hamiltonian += SS;
      std::cout << "Working.... " << i << "\n";
   }

   for (int i = 1; i <= L; ++i)
   {
      TotalSp = Sp(i);
      TotalSm = Sm(i);
      TotalSz = Sz(i);
   }

   TotalS2 = 0.5 * (prod(TotalSp, TotalSm, Ident) + prod(TotalSm, TotalSp, Ident))
      + prod(TotalSz, TotalSz, Ident);

   // Repeat for the auxiliary sites
   for (int i = 1; i < L; ++i)
   {
      MPOperator SS = 0.5 * Jperp * (prod(Sp(i, "aux"), Sm(i%L+1, "aux"), Ident)
                                     + prod(Sm(i, "aux"), Sp(i%L+1, "aux"), Ident))
         + Jz * prod(Sz(i), Sz(i%L+1), Ident);
      if (Beta != 0)
      {
         MPOperator SSIso = 0.5 * (prod(Sp(i, "aux"), Sm(i%L+1, "aux"), Ident)
                                   + prod(Sm(i, "aux"), Sp(i%L+1, "aux"), Ident))
            + prod(Sz(i, "aux"), Sz(i%L+1, "aux"), Ident);
         SS += Beta * prod(SSIso, SSIso, Ident);
      }

      AuxHamiltonian += SS;
      std::cout << "Working.... " << i << "\n";
   }

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);

   pheap::Initialize(argv[7], 1, 65536, 655360);
   pheap::ShutdownPersistent(OList);
}
