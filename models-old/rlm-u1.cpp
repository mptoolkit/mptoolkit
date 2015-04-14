// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "models/spinlessfermion-u1.h"

using std::string;

int main(int argc, char** argv)
{
   if (argc != 2)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: rlm-u1 <outfile> < couplings\n"
         //                << "L = number of lattice sites\n"
         //                << "t = hopping integral\n"
                << "outfile = file name for output lattice\n"
                << "couplings = epsilon_0, t_i epsilon_i for i = 1 .. L.\n";
      return 1;
   }

   std::string OutFile = argv[1];

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
   SiteBlock BathSite = CreateU1SpinlessFermion();
   SiteBlock ImpuritySite = CreateU1SpinlessFermion();

   // Create the lattice
   Lattice Bath = repeat(BathSite, L);
   Bath.fix_coordinates();
   Lattice Impurity = ImpuritySite;
   Impurity.fix_coordinates("0");
   Lattice MyLattice = join(Impurity, Bath);
   MyLattice.fix_coordinates_unique();

   // construct the operator list for the lattice
   OperatorList OpList(MyLattice);

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number

   // operators on the unwrapped lattice
   OperatorAtSite<OperatorList const, int> C(OpList, "C");
   OperatorAtSite<OperatorList const, int> CH(OpList, "CH");
   OperatorAtSite<OperatorList const, int> N(OpList, "N");
   MPOperator& Hamiltonian = OpList["H"];
   // hopping matrix elements
   for (int i = 0; i < L; ++i)
   {
      MPOperator Hopping 
         = prod(CH(i), C(i%L+1), Ident)
      	 - prod(C(i), CH(i%L+1), Ident);
      Hamiltonian += -Hop[i] * Hopping + Epsilon[i] * N(i);

      std::cout << "Working.... " << i << "\n";
   }
   // the final epsilon term
   Hamiltonian = Hamiltonian + Epsilon[L] * N(L);

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);

   pheap::Initialize(OutFile, 1, 65536, 655360);
   pheap::ShutdownPersistent(OList);
}
