// -*- C++ -*- $Id: siam-factorized-u1u1.cpp 514 2007-04-15 01:58:00Z ianmcc $

#include "pheap/pheap.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "models/spinlessfermion-u1.h"
#include "common/math_const.h"
#include <fstream>

using std::string;

int main(int argc, char** argv)
{
   if (argc != 7)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: siam-u1 <U> <Gamma12> <Gamma21> <bath1> <bath2> <outfile>\n"
                << "U = impurity U\n"
                << "Gamma12 = resonance width 1-2\n"
                << "Gamma21 = resonance width 2-1\n"
                << "bath1 = bath parameters for left chain\n"
                << "bath2 = bath parameters for right chain\n"
                << "outfile = file name for output lattice\n"
                      ;
      return 1;
   }

   double U = boost::lexical_cast<double>(argv[1]);
   double Gamma12 = boost::lexical_cast<double>(argv[2]);
   double Gamma21 = boost::lexical_cast<double>(argv[3]);
   std::string Bath1Name = argv[4];
   std::string Bath2Name = argv[5];
   std::string OutName = argv[6];

   TRACE(U)(Gamma12)(Gamma21);


   std::vector<double> Hop1, Epsilon1, Hop2, Epsilon2;
   double e,t;

   // Read Bath1
   std::ifstream In1(Bath1Name.c_str());
   In1 >> e;
   Epsilon1.push_back(e);
   while (In1 >> t >> e)
   {
      Hop1.push_back(t);
      Epsilon1.push_back(e);
   }

   // Read Bath2
   std::ifstream In2(Bath2Name.c_str());
   In2 >> e;
   Epsilon2.push_back(e);
   while (In2 >> t >> e)
   {
      Hop2.push_back(t);
      Epsilon2.push_back(e);
   }

   int L1 = Epsilon1.size();
   int L2 = Epsilon2.size();
   std::cout << "Lattice size = " << L1 << "+" << L2 << " = " << (L1+L2) << '\n';
   
   // Construct the site block
   SiteBlock Site = CreateU1SpinlessFermion();

   // Create the lattice
   Lattice LeftSite(Site);
   LeftSite.fix_coordinates("1");
   Lattice RightSite(Site);
   RightSite.fix_coordinates("2");

   Lattice LeftLattice = repeat(LeftSite, L1);
   LeftLattice.fix_coordinates_reverse();

   Lattice RightLattice = repeat(RightSite, L2);
   RightLattice.fix_coordinates();

   Lattice MyLattice = join(LeftLattice, RightLattice);
   MyLattice.fix_coordinates_unique();

   // construct the operator list for the lattice
   OperatorList OpList(MyLattice);

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number

   // operators on the unwrapped lattice
   OperatorAtSite<OperatorList const> C(OpList, "C");
   OperatorAtSite<OperatorList const> CH(OpList, "CH");
   OperatorAtSite<OperatorList const> N(OpList, "N");
   MPOperator& Hamiltonian = OpList["H"];

   // hopping matrix elements for bath 1
   for (int i = 1; i < L1; ++i)
   {
      Hamiltonian += -Hop1[i-1] * (CH(1,i) * C(1,i+1) - C(1,i) * CH(1,i+1))
         + Epsilon1[i-1] * N(1,i);
      std::cout << "Working.... " << i << "\n";
   }
   // the final epsilon term
   Hamiltonian += Epsilon1[L1-1] * N(1,L1);

   // hopping matrix elements for bath 2
   for (int i = 1; i < L2; ++i)
   {
      Hamiltonian += -Hop2[i-1] * (CH(2,i) * C(2,i+1) - C(2,i) * CH(2,i+1))
         + Epsilon2[i-1] * N(2,i);
      std::cout << "Working.... " << i << "\n";
   }
   // the final epsilon term
   Hamiltonian += Epsilon2[L2-1] * N(2,L2);

   // coulomb term
   Hamiltonian += U * N(1,1) * N(2,1);

   // Inter-bath couplings
   Hamiltonian -= std::sqrt(2 * Gamma12 / math_const::pi) * (CH(1,2)*C(2,1) - C(1,2)*CH(2,1));
   Hamiltonian -= std::sqrt(2 * Gamma21 / math_const::pi) * (CH(2,1)*C(1,2) - C(1,2)*CH(2,1));

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);
   pheap::Initialize(OutName, 1, 65536, 655360);
   pheap::ShutdownPersistent(OList);
}
