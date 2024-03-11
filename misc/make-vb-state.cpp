// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// misc/make-vb-state.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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
#include "interface/inittemp.h"

int main(int argc, char** argv)
{
   if (argc != 4)
   {
      std::cerr << "usage: make-vb-perm <lattice> <num-sites-in-omega> <output-operator-name>\n";
      return 1;
   }

   // load the lattice file
   pvalue_ptr<OperatorList> System = pheap::OpenPersistent(argv[1], mp_pheap::CacheSize());
   int NumSitesInOmega = boost::lexical_cast<int>(argv[2]);
   std::string OpName = argv[3];

   // get the lattice size
   int LatticeSize = System->size();

   // a shortcut to refer to the spin exchange operators
   OperatorAtSite<OperatorList const, int> Sp(*System, "Spp");
   OperatorAtSite<OperatorList const, int> Sm(*System, "Smm");

   // here, Left denotes the set in Omega, Right denotes \bar\Omega
   MPOperator SPlusLeft, SPlusRight, SMinusLeft, SMinusRight;

   for (int i = 1; i <= NumSitesInOmega; ++i)
   {
      SPlusLeft += Sp(i);
      SMinusLeft += Sm(i);
   }
   for (int i = NumSitesInOmega+1; i <= LatticeSize; ++i)
   {
      SPlusRight += Sp(i);
      SMinusRight += Sm(i);
   }

   // insert our operator into the lattice
   (*System.mutate())[OpName] = SPlusLeft*SMinusRight + SMinusLeft*SPlusRight;

   // save the lattice
   pheap::ShutdownPersistent(System);
}
