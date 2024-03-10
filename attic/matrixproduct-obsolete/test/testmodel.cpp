// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// attic/matrixproduct-obsolete/test/testmodel.cpp
//
// Copyright (C) 2004-2017 Ian McCulloch <ian@qusim.net>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER
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
