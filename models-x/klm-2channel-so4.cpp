// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "models/hubbard-so4.h"
#include "models/spin-su2.h"

typedef std::complex<double> complex;

int main(int argc, char** argv)
{
   if (argc != 6)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: kondo-2channel-so4 <L> <t> <J_K> <J_H> <J_Hboundary> <outfile>\n"
                << "L = number of lattice sites\n"
                << "t = hopping integral\n"
                << "J_K = Kondo coupling\n"
                << "J_H = Direct heisenberg term\n"
                << "J_Hboundary = Direct heisenberg term to the boundary spin 1/2\n"
                << "outfile = file name for output lattice.\n";
      return 1;
   }

   int L = boost::lexical_cast<int>(argv[1]);
   complex t = boost::lexical_cast<complex>(argv[2]);
   double Jk = boost::lexical_cast<double>(argv[3]);
   double Jh = boost::lexical_cast<double>(argv[4]);

   TRACE(L)(t)(Jk)(Jh);

   // Construct the site blocks
   SiteBlock bSiteA1 = CreateSO4HubbardSiteA("Q1", "S");
   SiteBlock bSiteB1 = CreateSO4HubbardSiteB("Q1", "S");
   SiteBlock bSiteA2 = CreateSO4HubbardSiteA("Q2", "S");
   SiteBlock bSiteB2 = CreateSO4HubbardSiteB("Q2", "S");
   SiteBlock bSiteS = CreateSU2SpinSite(1, "S");

   // construct the lattice
   Lattice UnitCellA(SymmetryList("Q1:SU(2),Q2:SU(2),S:SU(2)"), bSiteA1, bSiteS, bSiteA2);
   UnitCellA.fix_coordinates("1","0","2");
   Lattice UnitCellB(SymmetryList("Q1:SU(2),Q2:SU(2),S:SU(2)"), bSiteB1, bSiteS, bSiteB2);
   UnitCellB.fix_coordinates("1","0","2");

   Lattice MyLattice = repeat(Lattice(UnitCellA, UnitCellB), L/2);
   if (L%2 == 1)
      MyLattice = join(MyLattice, UnitCellA);
   MyLattice.fix_coordinates();

   // construct the operator list for the lattice
   OperatorList OpList(MyLattice);

   OperatorAtSite<OperatorList const, int, int> C(OpList, "C");
   OperatorAtSite<OperatorList const, int, int> CH(OpList, "CH");
   OperatorAtSite<OperatorList const, int, int> S(OpList, "S");
   OperatorAtSite<OperatorList const, int, int> N_S(OpList, "N_S");
   OperatorAtSite<OperatorList const, int, int> N_H(OpList, "N_H");
   OperatorAtSite<OperatorList const, int, int> Q1(OpList, "Q1");
   OperatorAtSite<OperatorList const, int, int> Q2(OpList, "Q2");
   MPOperator& Hamiltonian = OpList["H"];

   // Make some combined operators, for convenience
   OperatorAtSite<OperatorList, int> SiteS(OpList, "SiteS");
   OperatorAtSite<OperatorList, int> LocalS(OpList, "LocalS");
   OperatorAtSite<OperatorList, int> ConductionS(OpList, "ConductionS");
   OperatorAtSite<OperatorList, int> SiteN_S(OpList, "SiteN_S");
   OperatorAtSite<OperatorList, int> SiteN_H(OpList, "SiteN_H");
   OperatorAtSite<OperatorList, int> SiteQ1(OpList, "SiteQ1");
   OperatorAtSite<OperatorList, int> SiteQ2(OpList, "SiteQ2");

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number
   // hopping matrix elements

   complex HoppingValue = -2.0 * t;
   double SpinValue = -sqrt(3.0);

   for (int i = 1; i < L; ++i)
   {
      // Hopping in all bands
      for (int Band = 1; Band <= 2; ++Band)
      {
         MPOperator Hopping = HoppingValue * prod(C(Band,i), CH(Band,i%L+1), Ident);
         Hamiltonian += Hopping;
      }
      // Direct Heisenberg term
      if (Jh != 0) Hamiltonian += (Jh * SpinValue) * prod(S(0,i), S(0,i%L+1), Ident);

      std::cout << "Working.... " << i << "\n";
   }
   // Kondo coupling, and set the convenience operators
   for (int i = 1; i <= L; ++i)
   {
      ConductionS(i) = S(1,i) + S(2,i);

      Hamiltonian = Hamiltonian + (Jk * SpinValue) * prod(S(0,i), ConductionS(i), Ident);

      SiteS(i) = ConductionS(i) + S(0,i);
      LocalS(i) = S(0,i);
      SiteN_S(i) = N_S(1,i) + N_S(2,i);
      SiteN_H(i) = N_H(1,i) + N_H(2,i);
      SiteQ1(i) = Q1(1,i);
      SiteQ2(i) = Q2(2,i);

      std::cout << "Working.... " << i << "\n";
   }

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);

   pheap::Initialize(argv[5], 1, 65536, 655360);
   pheap::ShutdownPersistent(OList);
}
