// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "models/hubbard-so4.h"

typedef std::complex<double> complex;

int main(int argc, char** argv)
{
   if (argc != 5)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: hubbard-so4 <L> <t> <U> <outfile>\n"
                << "L = number of lattice sites\n"
                << "t = hopping integral\n"
                << "U = coupling constant\n"
                << "outfile = file name for output lattice.\n";
      return 1;
   }

   int L = boost::lexical_cast<int>(argv[1]);
   double t = boost::lexical_cast<double>(argv[2]);
   double U = boost::lexical_cast<double>(argv[3]);

   TRACE(L)(t)(U);

   Lattice UnitCell(CreateSO4HubbardSiteA(), CreateSO4HubbardSiteB());
   Lattice MyLattice = repeat(UnitCell, L/2);
   if (L%2 == 1)
      MyLattice = join(MyLattice, Lattice(CreateSO4HubbardSiteA()));

   MyLattice.fix_coordinates();

   // construct the operator list for the lattice
   OperatorList OpList(MyLattice);

   OperatorAtSite<OperatorList const, int> C(OpList, "C");
   OperatorAtSite<OperatorList const, int> CH(OpList, "CH");
   OperatorAtSite<OperatorList const, int> P(OpList, "P");
   OperatorAtSite<OperatorList, int> Bond(OpList, "Bond");
   MPOperator& Hamiltonian = OpList["H"];
   MPOperator& Hop = OpList["Hopping"];
   MPOperator& Coulomb = OpList["Coulomb"];
   MPOperator& Spin = OpList["Spin"];

   Hop = -1.0 * CreateRepeatedOperator(MyLattice, "CP", "CH");
   Coulomb = 0.25 * CreateRepeatedOperator(MyLattice, "P");
   Hamiltonian = t*Hop + U*Coulomb;
   Spin = -1.0 * CreateRepeatedOperator(MyLattice, "S", "S");

   // Bond operators
   for (int i = 1; i < L; ++i)
   {
      Bond(i) = -t * dot(C(i), CH(i%L+1)) + (U*0.125)*P(i) + (U*0.125)*P(i+1);
   }
   Bond(1) += (U*0.125)*P(1);
   Bond(L-1) += (U*0.125)*P(L);

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);

   pheap::Initialize(argv[4], 1, 65536, 655360);
   pheap::ShutdownPersistent(OList);
}
