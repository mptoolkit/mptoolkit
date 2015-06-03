// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "models/hubbard-u1u1.h"

int main(int argc, char** argv)
{
   if (argc != 4)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: mp-u1siam <U> <H> <outfile> < couplings\n"
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
      Epsilon.push_back(e);
      Hop.push_back(t);
   }

   int L = Hop.size();
   std::cout << "Lattice size (excluding impurity) = " << L << '\n';
   
   // Construct the site block
   SiteBlock Site = CreateU1HubbardSite();

   // construct a lattice of L copies of Site, for the bath
   Lattice Bath = repeat(Site, L);
   Bath.fix_coordinates();

   // The impurity site gets coordinate '0'
   Lattice Impurity(Site);
   Impurity.fix_coordinates("0");

   Lattice MyLattice = join(Impurity, Bath);
   MyLattice.fix_coordinates_unique();
   
   // construct the operator list for the lattice
   OperatorList OpList(MyLattice);

   OperatorAtSite<OperatorList const, int> CHup(OpList, "CHup");
   OperatorAtSite<OperatorList const, int> Cup(OpList, "Cup");
   OperatorAtSite<OperatorList const, int> CHdown(OpList, "CHdown");
   OperatorAtSite<OperatorList const, int> Cdown(OpList, "Cdown");
   OperatorAtSite<OperatorList const, int> N(OpList, "N");

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number
   MPOperator& Hamiltonian = OpList["H"];
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
   Hamiltonian += Epsilon[L] * N(L);

   // coulomb term
   Hamiltonian += U * OpList["Pdouble(0)"];

   // local field
   Hamiltonian += H0 * OpList["Sz(0)"];

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);

   pheap::Initialize(argv[3], 1, 65536, 655360);
   pheap::ShutdownPersistent(OList);
}
