// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "models/hubbard-u1su2.h"

typedef std::complex<double> complex;

int main(int argc, char** argv)
{
   if (argc != 3)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: hubbard-periodic-u1su2 <L> <outfile>\n"
                << "L = number of lattice sites\n"
                << "outfile = file name for output lattice.\n"
                << "\nOperators are:\n"
                << "H_t = symmetric (real) hopping (periodic)\n"
                << "H_j = antisymmetric (imaginary) hopping (periodic)\n"
                << "H_U = on-site Coulomb interaction (particle-hole symmetric)\n"
                << "H_J = nearest-neighbor spin exchange (periodic)\n"
                << "H_V = nearest-neighbor coulomb interaction (periodic)\n";
      return 1;
   }

   int L = boost::lexical_cast<int>(argv[1]);
   std::string OutFile = argv[2];

   // Construct the site block
   SiteBlock Site = CreateSU2HubbardSite();

   // construct a lattice of L copies of Site
   Lattice MyLattice = repeat(Site, L);
   MyLattice.fix_coordinates();

   // construct the operator list for the lattice
   OperatorList OpList(MyLattice);

   OperatorAtSite<OperatorList const, int> CH(OpList, "CH");
   OperatorAtSite<OperatorList const, int> C(OpList, "C");
   OperatorAtSite<OperatorList const, int> P(OpList, "P");
   OperatorAtSite<OperatorList const, int> N(OpList, "N");
   OperatorAtSite<OperatorList const, int> S(OpList, "S");
   OperatorAtSite<OperatorList const, int> LocalPg(OpList, "Pg");
   MPOperator& H_t = OpList["H_t"];
   MPOperator& H_j = OpList["H_j"];
   MPOperator& H_U = OpList["H_U"];
   MPOperator& H_J = OpList["H_J"];
   MPOperator& H_V = OpList["H_V"];

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number
   // hopping matrix elements

   // nearest neighbor terms
   for (int i = 1; i <= L; ++i)
   {
      H_t -= dot(CH(i), C(i%L+1)) + dot(C(i), CH(i%L+1));  // coefficient is -t
      H_j -= complex(0.0, 1.0) * (dot(CH(i), C(i%L+1)) - dot(C(i), CH(i%L+1)));  // coefficient is -t
      H_J -= dot(S(i), S(i%L+1));  // minus sign for SU(2) axial vector
      H_V += N(i)*N(i%L+1);
      std::cout << "Working.... " << i << "\n";
   }

   // on-site terms
   for (int i = 1; i <= L; ++i)
   {
      H_U += 0.25 * P(i);
      std::cout << "Working.... " << i << "\n";
   }

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);

   pheap::Initialize(OutFile, 1, 65536, 655360);
   pheap::ShutdownPersistent(OList);
}
