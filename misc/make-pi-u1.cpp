// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// misc/make-pi-u1.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Reseach publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "pheap/pheap.h"

int main(int argc, char** argv)
{
   if (argc != 2)
   {
      std::cerr << "usage: make-pi-u1 <lattice>\n";
      return 1;
   }

   // load the lattice file
   pvalue_ptr<OperatorList> System = pheap::OpenPersistent(argv[1], 655360);

   // get the lattice size
   int LatticeSize = System->size();

   // a shortcut to refer to the "S" (spin) operator
   OperatorAtSite<OperatorList const, int> Sp(*System, "Sp");
   OperatorAtSite<OperatorList const, int> Sm(*System, "Sm");
   OperatorAtSite<OperatorList const, int> Sz(*System, "Sz");

   // construct our output operator
   MPOperator Spip, Spim, Spiz;

   for (int i = 1; i <= LatticeSize; ++i)
   {
      Spip += pow(-1, i) * Sp(i);
      Spim += pow(-1, i) * Sm(i);
      Spiz += pow(-1, i) * Sz(i);
   }

   // insert our operator into the lattice
   (*System.mutate())["Spip"] = Spip;
   (*System.mutate())["Spim"] = Spim;
   (*System.mutate())["Spiz"] = Spiz;

   // save the lattice
   pheap::ShutdownPersistent(System);
}
