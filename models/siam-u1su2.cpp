// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "models/hubbard-u1su2.h"

int main(int argc, char** argv)
{
   if (argc != 3)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: siam-u1su2 <U> <outfile> < couplings\n"
                << "U = coupling constant\n"
                << "outfile = file name for output lattice\n"
                << "couplings = epsilon_0, t_i epsilon_i for i = 1 .. L.\n";
      return 1;
   }

   double U = boost::lexical_cast<double>(argv[1]);

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
   SiteBlock Site = CreateSU2HubbardSite();

   // construct a lattice of L copies of Site
   Lattice MyLattice = repeat(Site, L);
   MyLattice.fix_coordinates();

   // construct the operator list for the lattice
   OperatorList OpList(MyLattice);

   OperatorAtSite<OperatorList const, int> CH(OpList, "CH");
   OperatorAtSite<OperatorList const, int> C(OpList, "C");
   OperatorAtSite<OperatorList const, int> N(OpList, "N");

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number
   MPOperator& Hamiltonian = OpList["H"];
   // hopping matrix elements
   for (int i = 0; i < L; ++i)
   {
      MPOperator Hopping = prod(CH(i), C(i%L+1), Ident) + prod(C(i), CH(i%L+1), Ident);
      Hamiltonian += (-Hop[i] * std::sqrt(2.0)) * Hopping + Epsilon[i] * N(i);
   }
   // the final epsilon term
   Hamiltonian += Epsilon[L] * N(L);

   // coulomb term
   Hamiltonian += U * OpList["Pdouble(0)"];

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);

   pheap::Initialize(argv[2], 1, 65536, 655360);
   pheap::ShutdownPersistent(OList);
}
