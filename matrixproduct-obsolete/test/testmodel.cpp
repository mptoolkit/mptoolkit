// -*- C++ -*- $Id$
//
// A real test on the Hubbard model with SO(4) symmetry

#include "models/hubbard-so4.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"

void DoCheck(int L)
{
   double U = 5;
   Lattice UnitCell(CreateSO4HubbardSiteA(), CreateSO4HubbardSiteB());
   Lattice MyLattice = repeat(UnitCell, L/2);
   if (L%2 == 1)
      MyLattice = join(MyLattice, Lattice(CreateSO4HubbardSiteA()));

   OperatorList OpList(MyLattice);

   OperatorAtSite<OperatorList const, int> C(OpList, "C");
   OperatorAtSite<OperatorList const, int> CH(OpList, "CH");
   OperatorAtSite<OperatorList const, int> P(OpList, "P");

   // Construct a Hamiltonian the slow way
   MPOperator Hamiltonian;
   for (int i = 1; i < L; ++i)
   {
      Hamiltonian += -1.0 * dot(C(i), CH(i%L+1));
   }
   // coulomb repulsion
   for (int i = 1; i <= L; ++i)
   {
      Hamiltonian = Hamiltonian + (U/4.0) * P(i);
   }

   // Construct the same Hamiltonian as a periodic operator
   MPOperator Hamiltonian2 = (U/4.0) * CreatePeriodicOperator(MyLattice, "P")
      - 1.0*CreatePeriodicOperator(MyLattice, "CP", "CH");

   // Check that they are the same
   CHECK((Hamiltonian2-Hamiltonian).empty());
}

int main()
{
   for (int i = 1; i < 17; ++i)
   {
      DoCheck(i);
   }
}
