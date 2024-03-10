// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// misc/make-pi.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "pheap/pheap.h"

int main(int argc, char** argv)
{
   if (argc != 2)
   {
      std::cerr << "usage: make-pi <lattice>\n";
      return 1;
   }

   // load the lattice file
   pvalue_ptr<OperatorList> System = pheap::OpenPersistent(argv[1], 655360);

   // get the lattice size
   int LatticeSize = System->size();

   // a shortcut to refer to the "S" (spin) operator
   OperatorAtSite<OperatorList const, int> S(*System, "S");

   // construct our output operator
   MPOperator Spi;

   for (int i = 1; i <= LatticeSize; ++i)
   {
      Spi += pow(-1, i) * S(i);
   }

   // insert our operator into the lattice
   (*System.mutate())["Spi"] = Spi;

   // save the lattice
   pheap::ShutdownPersistent(System);
}
