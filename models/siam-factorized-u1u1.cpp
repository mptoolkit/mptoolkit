// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "models/factorized-u1u1.h"
#include "models/hubbard-u1u1.h"

using std::string;

int main(int argc, char** argv)
{
   if (argc != 4)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: siam-factorized-u1u1 <U> <H> <outfile> < couplings\n"
         //                << "L = number of lattice sites\n"
         //                << "t = hopping integral\n"
                << "U = impurity U\n"
                << "H = local magnetic field\n"
                << "outfile = file name for output lattice\n"
                << "couplings = epsilon_0, t_i epsilon_i for i = 1 .. L.\n";
      return 1;
   }

   double U = boost::lexical_cast<double>(argv[1]);
   double H0 = boost::lexical_cast<double>(argv[2]);

   std::vector<double> Epsilon;
   std::vector<double> Hop;

   double e,t;
   std::cin >> e;
   Epsilon.push_back(e);
   while (std::cin >> t >> e)
   {
      Hop.push_back(t);
      Epsilon.push_back(e);
   }

   int L = Hop.size();
   std::cout << "Lattice size (excluding impurity) = " << L << '\n';
   
   // Construct the site block
   SiteBlock DownSite = FactorizedSite(1, -0.5);
   SiteBlock UpSite = FactorizedSite(1, 0.5);
   SiteBlock ImpuritySite = CreateU1HubbardSite();

   // Create the lattice
   Lattice DownL(DownSite);
   DownL.fix_coordinates("down");
   Lattice UpL(UpSite);
   UpL.fix_coordinates("up");
   Lattice ImpurityL(ImpuritySite);
   ImpurityL.fix_coordinates("impurity");

   Lattice DownLattice = repeat(DownL, L);
   DownLattice.fix_coordinates_reverse();
   Lattice UpLattice = repeat(UpL, L);
   UpLattice.fix_coordinates();

   Lattice MyLattice = join(DownLattice, ImpurityL, UpLattice);
   MyLattice.fix_coordinates_unique();

   // construct the operator list for the lattice
   OperatorList OpList(MyLattice);

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number

   // operators on the unwrapped lattice
   OperatorAtSite<OperatorList const, string, int> C(OpList, "C");
   OperatorAtSite<OperatorList const, string, int> CH(OpList, "CH");
   OperatorAtSite<OperatorList const, string, int> LocalSz(OpList, "Sz");
   OperatorAtSite<OperatorList const, string, int> LocalN(OpList, "N");
   
   // make operators for the 'wrapped' lattice
   OperatorAtSite<OperatorList, int> CHup(OpList, "CHup");
   OperatorAtSite<OperatorList, int> Cup(OpList, "Cup");
   OperatorAtSite<OperatorList, int> CHdown(OpList, "CHdown");
   OperatorAtSite<OperatorList, int> Cdown(OpList, "Cdown");
   OperatorAtSite<OperatorList, int> N(OpList, "N");
   OperatorAtSite<OperatorList, int> Pdouble(OpList, "Pdouble");
   OperatorAtSite<OperatorList, int> Sz(OpList, "Sz");
   OperatorAtSite<OperatorList, int> Sp(OpList, "Sp");
   OperatorAtSite<OperatorList, int> Sm(OpList, "Sm");

   for (int i = 1; i <= L; ++i)
   {
      Cup(i) = C("up", i);
      CHup(i) = CH("up", i);
      Cdown(i) = C("down", i);
      CHdown(i) = CH("down", i);
      N(i) = LocalN("up", i) + LocalN("down", i);
      Pdouble(i) = prod(LocalN("up", i), LocalN("down", i), Ident);
      Sz(i) = LocalSz("up", i) + LocalSz("down", i);
      Sp(i) = prod(CH("up", i), C("down", i), 
                   QuantumNumber(MyLattice.GetSymmetryList(), "0,1"));
      Sm(i) = prod(CH("down", i), C("up", i),
                   QuantumNumber(MyLattice.GetSymmetryList(), "0,-1"));
      std::cout << "Working.... " << i << "\n";
   }
   // Site 0 is the impurity site
   Cup(0) = OpList["Cup(impurity)"];
   CHup(0) = OpList["CHup(impurity)"];
   Cdown(0) = OpList["Cdown(impurity)"];
   CHdown(0) = OpList["CHdown(impurity)"];
   N(0) = OpList["N(impurity)"];
   Pdouble(0) = OpList["Pdouble(impurity)"];
   Sz(0) = OpList["Sz(impurity)"];
   Sp(0) = OpList["Sp(impurity)"];
   Sm(0) = OpList["Sm(impurity)"];

   // Now we can refer to the operators as if they were on the original lattice

   MPOperator& Hamiltonian = OpList["H"];
   MPOperator& H_int = OpList["H_int"];
   
   // hopping matrix elements
   for (int i = 0; i < L; ++i)
   {
      MPOperator Hopping 
         = prod(CHup(i), Cup(i%L+1), Ident)
      	 - prod(Cup(i), CHup(i%L+1), Ident)
         + prod(CHdown(i), Cdown(i%L+1), Ident)
      	 - prod(Cdown(i), CHdown(i%L+1), Ident);
      Hamiltonian += -Hop[i] * Hopping + Epsilon[i] * N(i);

      std::cout << "Working.... " << i << "\n";
   }
   // the final epsilon term
   Hamiltonian = Hamiltonian + Epsilon[L] * N(L);

   // coulomb term
   Hamiltonian = Hamiltonian + U * OpList["Pdouble(0)"];
   H_int = U * OpList["Pdouble(0)"];

   // local field
   Hamiltonian = Hamiltonian - H0 * OpList["Sz(0)"];

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);

   pheap::Initialize(argv[3], 1, 65536, 655360);
   pheap::ShutdownPersistent(OList);
}
