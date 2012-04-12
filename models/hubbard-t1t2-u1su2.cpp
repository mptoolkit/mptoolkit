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
      std::cerr << "usage: hubbard-u1su2 <L> <outfile>\n"
                << "L = number of lattice sites\n"
                << "outfile = file name for output lattice.\n"
		<< "\nOutput operators are:\n"
		<< "H_t1 = nearest neighbor hopping\n"
		<< "H_t2 = next nearest neighbor hopping\n"
		<< "H_U = on-site coulomb repulsion (n_{up} * n_{down})\n"
		<< "H_U_sym = on-site coulomb repulsion (particle-hole symmetric alternative)\n"
		<< "H_V = nearest neighbor coulomb repulsion\n"
		<< "H_J1 = nearest neighbor spin-spin interaction\n"
		<< "H_J2 = nearest neighbor spin-spin interaction\n"
	 ;
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
   MPOperator& Hop1 = OpList["H_t1"];
   MPOperator& Hop2 = OpList["H_t2"];
   MPOperator& CoulombU_sym = OpList["H_U_sym"];
   MPOperator& CoulombU = OpList["H_U"];
   MPOperator& CoulombV = OpList["H_V"];
   MPOperator& Spin1 = OpList["H_J1"];
   MPOperator& Spin2 = OpList["H_J2"];

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number
   // hopping matrix elements

   double Sqrt2 = -std::sqrt(2.0);  // negative because the term is -t * ....
   double Sqrt3 = -std::sqrt(3.0);  // negative because of SU(2) vector phase

   Hop1 = CreateRepeatedOperator(MyLattice, "CP", "CH")
      + CreateRepeatedOperator(MyLattice, "CHP", "C"); // no need for - sign, already -ve

   Hop2 = CreateRepeatedOperator(MyLattice, "CP", "P", "CH") 
      + CreateRepeatedOperator(MyLattice, "CHP", "P", "C");

   CoulombU_sym = 0.25 * CreateRepeatedOperator(MyLattice, "P");
   CoulombU = CreateRepeatedOperator(MyLattice, "Pdouble");
   
   CoulombV = CreateRepeatedOperator(MyLattice, "N", "N");

   Spin1 = -1.0 * CreateRepeatedOperator(MyLattice, "S", "S");
   Spin2 = -1.0 * CreateRepeatedOperator(MyLattice, "S", "I", "S");

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);

   pheap::Initialize(OutFile, 1, 65536, 655360);
   pheap::ShutdownPersistent(OList);
}
