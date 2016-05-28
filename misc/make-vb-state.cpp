// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// misc/make-vb-state.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
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
