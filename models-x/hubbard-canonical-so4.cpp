// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "models/hubbard-so4.h"

typedef std::complex<double> complex;

int main(int argc, char** argv)
{
   if (argc != 5)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: hubbard-canonical-so4 <L> <t> <U> <outfile>\n"
                << "L = number of lattice sites\n"
                << "t = hopping integral\n"
                << "U = coupling constant\n"
                << "outfile = file name for output lattice.\n";
      return 1;
   }

   int L = boost::lexical_cast<int>(argv[1]);
   complex t = boost::lexical_cast<complex>(argv[2]);
   double U = boost::lexical_cast<double>(argv[3]);

   TRACE(L)(t)(U);

   // Construct the site blocks
   SiteBlock RealSiteA = CreateSO4HubbardSiteA();
   SiteBlock RealSiteB = CreateSO4HubbardSiteB();
   SiteBlock AuxSiteA = CreateSO4HubbardSiteA("Qaux", "Saux");
   SiteBlock AuxSiteB = CreateSO4HubbardSiteB("Qaux", "Saux");

   // construct the lattice
   Lattice UnitCellA(SymmetryList("Q:SU(2),S:SU(2),Qaux:SU(2),Saux:SU(2)"), RealSiteA, AuxSiteA);
   UnitCellA.fix_coordinates("", "aux");
   Lattice UnitCellB(SymmetryList("Q:SU(2),S:SU(2),Qaux:SU(2),Saux:SU(2)"), RealSiteB, AuxSiteB);
   UnitCellB.fix_coordinates("", "aux");
   Lattice MyLattice = repeat(join(UnitCellA, UnitCellB), L/2);
   if (L%2 == 1)
      MyLattice = join(MyLattice, UnitCellA);
   MyLattice.fix_coordinates_prepend();

   // construct the operator list for the lattice
   OperatorList OpList(MyLattice);

   OperatorAtSite<OperatorList const> C(OpList, "C");
   OperatorAtSite<OperatorList const> CH(OpList, "CH");
   OperatorAtSite<OperatorList const> P(OpList, "P");
   MPOperator& Hamiltonian = OpList["H"];
   MPOperator& AuxHamiltonian = OpList["AuxH"];
   MPOperator& SplitHamiltonian = OpList["SplitH"];
   MPOperator& Hop = OpList["Hopping"];
   MPOperator& AuxHop = OpList["AuxHopping"];
   MPOperator& SplitHop = OpList["SplitHopping"];
   MPOperator& Coulomb = OpList["Coulomb"];
   MPOperator& AuxCoulomb = OpList["AuxCoulomb"];
   MPOperator& SplitCoulomb = OpList["SplitCoulomb"];

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number

   // hopping matrix elements
   for (int i = 1; i < L; ++i)
   {
      MPOperator Hopping = -2.0 * prod(C(i), CH(i%L+1), Ident);
      Hop += Hopping;
      Hopping *= t;
      Hamiltonian += Hopping;
      std::cout << "Working.... " << i << "\n";
   }
   // coulomb repulsion
   for (int i = 1; i <= L; ++i)
   {
      Hamiltonian = Hamiltonian + (U/4.0) * P(i);
      Coulomb += 0.25 * P(i);
   }
   // repeat the same for the auxiliary sites
   // hopping matrix elements
   for (int i = 1; i < L; ++i)
   {
      MPOperator Hopping = -2.0 * prod(C(i,"aux"), CH(i%L+1,"aux"), Ident);
      AuxHop += Hopping;
      Hopping *= t;
      AuxHamiltonian += Hopping;
      std::cout << "Working.... " << i << "\n";
   }
   // coulomb repulsion
   for (int i = 1; i <= L; ++i)
   {
      AuxHamiltonian += (U/4.0) * P(i,"aux");
      AuxCoulomb += 0.25 * P(i);
      std::cout << "Working.... " << i << "\n";
   }   

   // Split versions
   SplitHamiltonian = 0.5 * (Hamiltonian - AuxHamiltonian);
   SplitCoulomb = 0.5 * (Coulomb - AuxCoulomb);
   SplitHop = 0.5 * (Hop - AuxHop);

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);

   pheap::Initialize(argv[4], 1, 65536, 655360);
   pheap::ShutdownPersistent(OList);
}
