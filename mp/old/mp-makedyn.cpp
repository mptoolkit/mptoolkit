// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-makedyn.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include <iostream>
#include "mp/copyright.h"

int main(int argc, char** argv)
{
   if (argc != 5)
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "usage: mp-makedyn <lattice> <groundstate-energy> <frequency> <broadening>\n";
      return 1;
   }

   pvalue_ptr<OperatorList> System = pheap::OpenPersistent(argv[1], 655360);
   TRACE("here");
   //   pvalue_ptr<MPWavefunction> Psi = pheap::ImportHeap(argv[2]);
   double Freq = boost::lexical_cast<double>(argv[3]);
   double Broad = boost::lexical_cast<double>(argv[4]);

   MPOperator Ham = (*System)["H"];

   MPOperator Ident = (*System)["I"];

   double Energy = boost::lexical_cast<double>(argv[2]); // expectation(*Psi, Ham, *Psi).real();

   std::string OutName = std::string("G_H(") + argv[3] + "," + argv[4] + ")";
   (*System.mutate())[OutName] = std::complex<double>(Energy + Freq, Broad) * Ident - Ham;

   OutName = std::string("G2_H(") + argv[3] + "," + argv[4] + ")";
   MPOperator Part = (Energy + Freq) * Ident - Ham;
   (*System.mutate())[OutName] = prod(Part, Part, Part.TransformsAs()) + (Broad*Broad) * Ident;

   pheap::ShutdownPersistent(System);
}
