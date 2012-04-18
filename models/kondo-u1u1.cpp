// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "models/hubbard-u1u1.h"
#include "models/spin-u1.h"

int main(int argc, char** argv)
{
   if (argc != 5)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: kondo-u1u1 <L> <first-impurity-location> <second-impurity-location> <outfile>\n"
                << "L = number of lattice sites\n"
		<< "impurity-location = location of the impurity\n"
                << "outfile = file name for output lattice.\n"
		<< "\nOperators are:\n"
		<< "H_t = nearest neightbor hopping\n"
		<< "H_J1 = Kondo interaction to the first impurity (+ve J is antiferromagnetic coupling) \n"
		<< "H_J2 = Kondo interaction to the second impurity\n"
		<< "J_J1J2 = Direct spin-spin interaction between impurities\n"
	 ;
      return 1;
   }

   int L = boost::lexical_cast<int>(argv[1]);
   int Loc1 = boost::lexical_cast<int>(argv[2]);
   int Loc2 = boost::lexical_cast<int>(argv[3]);
   std::string OutFile = argv[4];
   
   // Construct the site blocks
   SiteBlock Site = CreateU1HubbardSite();
   SiteBlock Impurity = CreateU1SpinSite(0.5);

   Lattice LeftPart = repeat(Site, Loc1);
   LeftPart.fix_coordinates_prepend();

   Lattice Imp1 = Lattice("N:U(1),Sz:U(1)", Impurity);
   Imp1.fix_coordinates("s1");
   
   Lattice CentralPart = repeat(Site, Loc2-Loc1);
   CentralPart.fix_coordinates_prepend_starting_from(Loc1+1);

   Lattice Imp2 = Lattice("N:U(1),Sz:U(1)", Impurity);
   Imp2.fix_coordinates("s2");

   Lattice RightPart = repeat(Site, L - Loc2);
   RightPart.fix_coordinates_prepend_starting_from(Loc2+1);

   Lattice MyLattice = join(LeftPart, Imp1, CentralPart, Imp2, RightPart);
   MyLattice.fix_coordinates_unique();

   // for (int i = 1; i <= MyLattice.size(); ++i)
   //{
   //    std::cout << MyLattice.coordinate_at_site(i) << '\n';
   // }

   // construct the operator list for the lattice
   OperatorList OpList(MyLattice);

   OperatorAtSite<OperatorList const> CHup(OpList, "CHup");
   OperatorAtSite<OperatorList const> Cup(OpList, "Cup");
   OperatorAtSite<OperatorList const> CHdown(OpList, "CHdown");
   OperatorAtSite<OperatorList const> Cdown(OpList, "Cdown");
   OperatorAtSite<OperatorList const> P(OpList, "P");
   OperatorAtSite<OperatorList const> Sz(OpList, "Sz");
   OperatorAtSite<OperatorList const> Sp(OpList, "Sp");
   OperatorAtSite<OperatorList const> Sm(OpList, "Sm");
   MPOperator& H_t = OpList["H_t"];
   MPOperator& H_J1 = OpList["H_J1"];
   MPOperator& H_J2 = OpList["H_J2"];
   MPOperator& H_J1J2 = OpList["H_J1J2"];

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number
   // hopping matrix elements
   for (int i = 1; i < L; ++i)
   {
      H_t += -1.0 * (
		     prod(CHup(i), Cup(i%L+1), Ident)
		     - prod(Cup(i), CHup(i%L+1), Ident)
		     + prod(CHdown(i), Cdown(i%L+1), Ident)
		     - prod(Cdown(i), CHdown(i%L+1), Ident)
		     );
      std::cout << "Working.... " << i << "\n";
   }

   // Kondo interactions
   H_J1 = Sz(Loc1)*Sz("s1") + 0.5*(Sp(Loc1)*Sm("s1") + Sm(Loc1)*Sp("s1"));
   H_J2 = Sz(Loc2)*Sz("s2") + 0.5*(Sp(Loc2)*Sm("s2") + Sm(Loc2)*Sp("s2"));
   H_J1J2 = Sz("s1")*Sz("s2") + 0.5*(Sp("s1")*Sm("s2") + Sm("s1")*Sp("s2"));

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);

   pheap::Initialize(OutFile, 1, 65536, 655360);
   pheap::ShutdownPersistent(OList);
}
