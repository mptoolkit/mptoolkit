// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "models/hubbard-u1u1.h"

int main(int argc, char** argv)
{
   if (argc != 5)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: hubbard-canonical-u1u1 <L> <t> <U> <outfile>\n"
                << "L = number of lattice sites\n"
                << "t = hopping integral\n"
                << "U = coupling constant\n"
                << "outfile = file name for output lattice.\n";
      return 1;
   }

   int L = boost::lexical_cast<int>(argv[1]);
   double t = boost::lexical_cast<double>(argv[2]);
   double U = boost::lexical_cast<double>(argv[3]);

   TRACE(L)(t)(U);
   
   // Construct the site block
   SiteBlock RealSite = CreateU1HubbardSite();
   SiteBlock AuxSite = CreateU1HubbardSite("Naux", "Szaux");

   // construct a lattice of L copies of Site
   Lattice UnitCell(SymmetryList("N:U(1),Sz:U(1),Naux:U(1),Szaux:U(1)"), RealSite, AuxSite);
   UnitCell.fix_coordinates("", "aux");
   Lattice MyLattice = repeat(UnitCell, L);
   MyLattice.fix_coordinates_prepend();

   // construct the operator list for the lattice
   OperatorList OpList(MyLattice);

   OperatorAtSite<OperatorList const> CHup(OpList, "CHup");
   OperatorAtSite<OperatorList const> Cup(OpList, "Cup");
   OperatorAtSite<OperatorList const> CHdown(OpList, "CHdown");
   OperatorAtSite<OperatorList const> Cdown(OpList, "Cdown");
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
      MPOperator Hopping = -1.0 * (
				     prod(CHup(i), Cup(i%L+1), Ident)
				   - prod(Cup(i), CHup(i%L+1), Ident)
				   + prod(CHdown(i), Cdown(i%L+1), Ident)
				   - prod(Cdown(i), CHdown(i%L+1), Ident)
				 );
      Hop += Hopping;
      Hopping *= t;
      Hamiltonian += Hopping;
      std::cout << "Working.... " << i << "\n";
   }
   // coulomb repulsion
   for (int i = 1; i <= L; ++i)
   {
      Coulomb += 0.25 * P(i);
      Hamiltonian += (U/4.0) * P(i);
      std::cout << "Working.... " << i << "\n";
   }

   // repeat the same for the auxiliary sites
   // hopping matrix elements
   std::cout << "Auxiliary Hamiltonian...\n";
   for (int i = 1; i < L; ++i)
   {
      MPOperator Hopping = -1.0 * (
				     prod(CHup(i,"aux"), Cup(i%L+1,"aux"), Ident)
				   - prod(Cup(i,"aux"), CHup(i%L+1, "aux"), Ident)
				   + prod(CHdown(i,"aux"), Cdown(i%L+1,"aux"), Ident)
				   - prod(Cdown(i,"aux"), CHdown(i%L+1,"aux"), Ident)
				 );
      AuxHop += Hopping;
      Hopping *= t;
      AuxHamiltonian += Hopping;
      std::cout << "Working.... " << i << "\n";
   }
   // coulomb repulsion
   for (int i = 1; i <= L; ++i)
   {
      AuxCoulomb += 0.25 * P(i,"aux");
      AuxHamiltonian += (U/4.0) * P(i,"aux");
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
