// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "models/hubbard-u1u1.h"

int main(int argc, char** argv)
{
   if (argc != 3)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: hubbard-u1u1 <L> <outfile>\n"
                << "L = number of lattice sites\n"
                << "outfile = file name for output lattice.\n"
		<< "\nOperators are:\n"
		<< "H_tup = hopping for the up spins\n"
		<< "H_tdown = hopping for the down spins\n"
                << "H_t = spin-symmetric hopping (equal to H_tup + H_tdown)\n"
                << "H_U = on-site Coulomb interaction (particle-hole symmetric)\n"
                << "H_V = nearest-neighbor coulomb interaction\n";
      return 1;
   }

   int L = boost::lexical_cast<int>(argv[1]);
   std::string FName = argv[2];

   // Construct the site block
   SiteBlock Site = CreateU1HubbardSite();

   // construct a lattice of L copies of Site
   Lattice MyLattice(L, Site);
   MyLattice.fix_coordinates();

   // construct the operator list for the lattice
   OperatorList OpList(MyLattice);

   OperatorAtSite<OperatorList const, int> CHup(OpList, "CHup");
   OperatorAtSite<OperatorList const, int> Cup(OpList, "Cup");
   OperatorAtSite<OperatorList const, int> CHdown(OpList, "CHdown");
   OperatorAtSite<OperatorList const, int> Cdown(OpList, "Cdown");
   OperatorAtSite<OperatorList const, int> P(OpList, "P");
   OperatorAtSite<OperatorList const, int> R(OpList, "R");
   OperatorAtSite<OperatorList const, int> N(OpList, "N");
   OperatorAtSite<OperatorList const, int> Sp(OpList, "Sp");
   OperatorAtSite<OperatorList const, int> Sm(OpList, "Sm");
   OperatorAtSite<OperatorList const, int> Sz(OpList, "Sz");
   MPOperator& H_tup = OpList["H_tup"];
   MPOperator& H_tdown = OpList["H_tdown"];
   MPOperator& H_t = OpList["H_t"];
   MPOperator& H_U = OpList["H_U"];
   MPOperator& H_V = OpList["H_V"];

   MPOperator& TotalR = OpList["R"]; // spatial reflection, ie \up\down -> -\up\down
   MPOperator& R_A = OpList["R_A"];
   MPOperator& R_B = OpList["R_B"];

   MPOperator& TotalSp = OpList["Sp"]; 
   MPOperator& TotalSm = OpList["Sm"];
   MPOperator& TotalSz = OpList["Sz"];

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number
   // hopping matrix elements
   for (int i = 1; i < L; ++i)
   {
      H_tup -= dot(CHup(i), Cup(i%L+1)) - dot(Cup(i), CHup(i%L+1));
      H_tdown -= dot(CHdown(i), Cdown(i%L+1)) - dot(Cdown(i), CHdown(i%L+1));
      H_V += N(i)*N(i%L+1);
      std::cout << "Working.... " << i << "\n";
   }
   // coulomb repulsion and sublattice reflection
   TotalR = R_A = R_B = OpList["I"];  // multiplicative identity
   for (int i = 1; i <= L; ++i)
   {
      H_U += 0.25 * P(i);
      TotalR = TotalR * R(i);
      if (i % 2 == 1)
	 R_A = R_A * R(i);
      else
	 R_B = R_B * R(i);

      TotalSp += Sp(i);
      TotalSm += Sm(i);
      TotalSz += Sz(i);

      std::cout << "Working.... " << i << "\n";
   }
   H_t = H_tup + H_tdown;

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);

   pheap::Initialize(FName, 1, 65536, 655360);
   pheap::ShutdownPersistent(OList);
}
