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
      std::cerr << "usage: klm-so4 <L> <t> <J_K> <J_H> <outfile>\n"
                << "L = number of lattice sites\n"
                << "t = hopping integral\n"
                << "J_K = Kondo coupling\n"
                << "J_H = Direct heisenberg term\n"
                << "outfile = file name for output lattice.\n";
      return 1;
   }

   int L = boost::lexical_cast<int>(argv[1]);
   complex t = boost::lexical_cast<complex>(argv[2]);
   double Jk = boost::lexical_cast<double>(argv[3]);
   double Jh = boost::lexical_cast<double>(argv[4]);

   TRACE(L)(t)(Jk)(Jh);

   // Construct the site blocks
   SiteBlock bSiteA = CreateSO4HubbardSiteA("Q", "S");
   SiteBlock bSiteB = CreateSO4HubbardSiteB("Q", "S");
   SiteBlock bSiteS = CreateSU2SpinSite(0.5, "S");

   // construct the lattice
   Lattice UnitCell(SymmetryList("Q:SU(2),S:SU(2)"), bSiteS, bSiteA, bSiteS, bSiteB);
   UnitCell.fix_coordinates("f","c","f","c");

   Lattice MyLattice = repeat(UnitCell, L/2);
   if (L%2 == 1)
   {
      Lattice HalfCell(SymmetryList("Q:SU(2),S:SU(2)"), bSiteS, bSiteA);
      HalfCell.fix_coordinates("f", "c");
      MyLattice = join(MyLattice, HalfCell);
   }

   MyLattice.fix_coordinates();

   // construct the operator list for the lattice
   OperatorList OpList(MyLattice);

   OperatorAtSite<OperatorList> C(OpList, "C");
   OperatorAtSite<OperatorList> CH(OpList, "CH");
   OperatorAtSite<OperatorList const> S(OpList, "S");
   OperatorAtSite<OperatorList const> N_S(OpList, "N_S");
   OperatorAtSite<OperatorList const> N_H(OpList, "N_H");
   OperatorAtSite<OperatorList const> Q(OpList, "Q");
   MPOperator& Hamiltonian = OpList["H"];
   OperatorAtSite<OperatorList, int> fS(OpList, "fS");
   OperatorAtSite<OperatorList, int> cS(OpList, "cS");

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number
   // hopping matrix elements

   complex HoppingValue = -2.0 * t;
   double SpinValue = -sqrt(3.0);

   // operators acting as if we had a conduction and spin rolled into a supersite
   for (int i = 1; i <= L; ++i)
   {
      fS(i) = S("f", i);
      cS(i) = S("c", i);
      C(i) = C("c", i);
      CH(i) = CH("c", i);
   }

   for (int i = 1; i < L; ++i)
   {
      // Hopping
      MPOperator Hopping = HoppingValue * prod(C("c",i), CH("c",i%L+1), Ident);
      Hamiltonian += Hopping;
      
      // Direct Heisenberg term
      if (Jh != 0) Hamiltonian += (Jh * SpinValue) * prod(S("f",i), S("f",i%L+1), Ident);

      std::cout << "Working.... " << i << "\n";
   }
   // Kondo coupling, and set the convenience operators
   for (int i = 1; i <= L; ++i)
   {
      Hamiltonian = Hamiltonian + (Jk * SpinValue) * prod(S("f",i), S("c",i), Ident);
      
      std::cout << "Working.... " << i << "\n";
   }

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);

   pheap::Initialize(argv[5], 1, 65536, 655360);
   pheap::ShutdownPersistent(OList);
}
