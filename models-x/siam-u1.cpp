// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "models/spinlessfermion-u1.h"
//#include "models/hardcoreboson-u1.h"
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
   MPOperator& H1 = OpList["H1"];
   MPOperator& H2 = OpList["H2"];
   MPOperator& HU = OpList["HU"];
   MPOperator& HGamma12 = OpList["HGamma12"];
   MPOperator& HGamma21 = OpList["HGamma21"];
   MPOperator& HGamma12Bar = OpList["HGamma12Bar"];
   MPOperator& HGamma21Bar = OpList["HGamma21Bar"];
   MPOperator& Hamiltonian = OpList["H"];

   MPOperator& H1Bar = OpList["H1Bar"];
   MPOperator& H2Bar = OpList["H2Bar"];

   // hopping matrix elements for bath 1
   for (int i = 1; i < L1; ++i)
   {
      H1 += -Hop1[i-1] * (CH(1,i)*C(1,i+1) + CH(1,i+1)*C(1,i))
         + Epsilon1[i-1] * N(1,i);

      H1Bar += -Hop1[i-1] * (CH(1,i) * C(1,i+1) - C(1,i) * CH(1,i+1))
         + Epsilon1[i-1] * N(1,i);

      std::cout << "Working.... " << i << "\n";
   }
   // the final epsilon term
   H1 += Epsilon1[L1-1] * N(1,L1);

   H1Bar += Epsilon1[L1-1] * N(1,L1);

   // hopping matrix elements for bath 2
   for (int i = 1; i < L2; ++i)
   {
      H2 += -Hop2[i-1] * (CH(2,i)*C(2,i+1) + CH(2,i+1)*C(2,i))
         + Epsilon2[i-1] * N(2,i);

      H2Bar += -Hop2[i-1] * (CH(2,i) * C(2,i+1) - C(2,i) * CH(2,i+1))
         + Epsilon2[i-1] * N(2,i);

      std::cout << "Working.... " << i << "\n";
   }
   // the final epsilon term
   H2 += Epsilon2[L2-1] * N(2,L2);

   H2Bar += Epsilon2[L2-1] * N(2,L2);

   // coulomb term
   HU = U * N(1,1) * N(2,1);

   // Inter-bath couplings
   HGamma12 = -std::sqrt(2.0 * Gamma12 / math_const::pi) * (CH(1,2)*C(2,1) + CH(2,1)*C(1,2));
   HGamma21 = -std::sqrt(2.0 * Gamma21 / math_const::pi) * (CH(1,1)*C(2,2) + CH(2,2)*C(1,1));  

   HGamma12Bar = -std::sqrt(2.0 * Gamma12 / math_const::pi) * (CH(1,2)*C(2,1) - C(1,2)*CH(2,1));
   HGamma21Bar = -std::sqrt(2.0 * Gamma21 / math_const::pi) * (CH(1,1)*C(2,2) - C(1,1)*CH(2,2));  

   Hamiltonian = H1 + H2 + HU + HGamma12 + HGamma21;











   MPOperator& HWeird = OpList["Hweird"];

   // hopping matrix elements for bath 1
   for (int i = 1; i < L1; ++i)
   {
      HWeird += -Hop1[i-1] * (CH(1,i) * C(1,i+1) - C(1,i) * CH(1,i+1))
         + Epsilon1[i-1] * N(1,i);
      std::cout << "Working.... " << i << "\n";
   }
   // the final epsilon term
   HWeird += Epsilon1[L1-1] * N(1,L1);

   // hopping matrix elements for bath 2
   for (int i = 1; i < L2; ++i)
   {
      HWeird += -Hop2[i-1] * (CH(2,i) * C(2,i+1) - C(2,i) * CH(2,i+1))
         + Epsilon2[i-1] * N(2,i);
      std::cout << "Working.... " << i << "\n";
   }
   // the final epsilon term
   HWeird += Epsilon2[L2-1] * N(2,L2);

   // coulomb term
   HWeird += U * N(1,1) * N(2,1);

   // Inter-bath couplings
   HWeird -= std::sqrt(2 * Gamma12 / math_const::pi) * (CH(1,2)*C(2,1) - C(1,2)*CH(2,1));
   HWeird -= std::sqrt(2 * Gamma21 / math_const::pi) * (CH(1,1)*C(2,2) - C(1,1)*CH(2,2));



















   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);
   pheap::Initialize(OutName, 1, 65536, 655360);
   pheap::ShutdownPersistent(OList);
}
