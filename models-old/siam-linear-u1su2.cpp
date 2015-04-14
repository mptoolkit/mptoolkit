// -*- C++ -*- $Id$

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "models/hubbard-u1su2.h"
#include "mp/copyright.h"
#include "pstream/pfilestream.h"

int main(int argc, char** argv)
{
   if (argc != 6)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: mp-su2siam <U> <V> <e0> <L> <outfile>\n"
         //                << "L = number of lattice sites\n"
         //                << "t = hopping integral\n"
                << "U = coupling constant\n"
		<< "V = hybridization\n"
		<< "e0 = impurity level"
		<< "L = lattice size (excluding impurity)"
                << "outfile = file name for output lattice\n";
      return 1;
   }

   double U = boost::lexical_cast<double>(argv[1]);
   double V = boost::lexical_cast<double>(argv[2]);
   double e0 = boost::lexical_cast<double>(argv[3]);
   int L = boost::lexical_cast<int>(argv[4]);

   std::vector<double> Epsilon;
   std::vector<double> Hop;

   Epsilon.push_back(e0);
   Epsilon.push_back(0.0);
   Hop.push_back(1.0);
   for (int i = 1; i < L; ++i)
   {
      Epsilon.push_back(0.0);
      Hop.push_back(1.0);
   }

   // Construct the site block
   SiteBlock Site = CreateSU2HubbardSite();

   // construct a lattice of L copies of Site
   Lattice MyLattice;
   for (int i = 0; i <= L; ++i)
   {
      MyLattice.Append(boost::lexical_cast<std::string>(i), Site);
   }

   // construct the operator list for the lattice
   OperatorList OpList(MyLattice);

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number
   MPOperator Hamiltonian;
   MPOperator E, F;
   // hopping matrix elements
   for (int i = 0; i < L; ++i)
   {
      MPOperator Hopping 
         = prod(OpList["CH("+boost::lexical_cast<std::string>(i)+")"],
                OpList["C("+boost::lexical_cast<std::string>(i%L+1)+")"],
                Ident)
      	 + prod(OpList["C("+boost::lexical_cast<std::string>(i)+")"],
      		OpList["CH("+boost::lexical_cast<std::string>(i%L+1)+")"],
		Ident);
      Hamiltonian = Hamiltonian + (Hop[i] * std::sqrt(2.0)) * Hopping;
      Hamiltonian = Hamiltonian + Epsilon[i] * OpList["N("+boost::lexical_cast<std::string>(i)+")"];

      std::cout << "Working.... " << i << std::endl;
   }
   // the final epsilon term
   Hamiltonian = Hamiltonian + Epsilon[L] * OpList["N("+boost::lexical_cast<std::string>(L)+")"];

   // coulomb term
   Hamiltonian = Hamiltonian + U * OpList["Pdouble(0)"];

   // save the Hamiltonian into the operator list
   OpList["H"] = Hamiltonian;

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);

   pheap::Initialize(argv[5], 1, 65536, 655360);
   pheap::ShutdownPersistent(OList);
}
