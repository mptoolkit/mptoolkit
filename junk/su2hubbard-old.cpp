// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// junk/su2hubbard-old.cpp
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
#include "models/su2hubbard.h"
#include "matrixproduct/copyright.h"
#include "pstream/pfilestream.h"

int main(int argc, char** argv)
{
   if (argc != 6)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: mp-su2hubbard <L> <t-real> <t-imag> <U> <outfile>\n"
                << "L = number of lattice sites\n"
                << "t = hopping integral\n"
                << "U = coupling constant\n"
                << "outfile = file name for output lattice.\n";
      return 1;
   }

   int L = boost::lexical_cast<int>(argv[1]);
   std::complex<double> t(boost::lexical_cast<double>(argv[2]), boost::lexical_cast<double>(argv[3]));
   double U = boost::lexical_cast<double>(argv[4]);

   TRACE(L)(t)(U);

   // Construct the site block
   SiteBlock Site = CreateSU2HubbardSite();

   // construct a lattice of L copies of Site
   Lattice MyLattice;
   for (int i = 1; i <= L; ++i)
   {
      MyLattice.Append(boost::lexical_cast<std::string>(i), Site);
   }

   // construct the operator list for the lattice
   OperatorList OpList(MyLattice);

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number
   MPOperator Hamiltonian;
   MPOperator E, F;
   // hopping matrix elements
   for (int i = 1; i < L; ++i)
   {
      MPOperator Hopping
         = (std::sqrt(2.0) * t) * prod(OpList["CH("+boost::lexical_cast<std::string>(i)+")"],
                OpList["C("+boost::lexical_cast<std::string>(i%L+1)+")"],
                Ident)
         + (std::sqrt(2.0) * conj(t)) * prod(OpList["C("+boost::lexical_cast<std::string>(i)+")"],
                OpList["CH("+boost::lexical_cast<std::string>(i%L+1)+")"],
                Ident);
      Hamiltonian = Hamiltonian + Hopping;

      if (i%2 == 1)
         E = E + Hopping;
      else
         F = F + Hopping;

      std::cout << "Working.... " << i << "\n";
   }
   // coulomb repulsion
   for (int i = 1; i <= L; ++i)
   {
      Hamiltonian = Hamiltonian + (U/4.0) * OpList["P("+boost::lexical_cast<std::string>(i)+")"];
      //      E = E + (U/8.0) * OpList["P("+boost::lexical_cast<std::string>(i)+")"];
      //      F = F + (U/8.0) * OpList["P("+boost::lexical_cast<std::string>(i)+")"];
      std::cout << "Working.... " << i << "\n";
   }

   // save the Hamiltonian into the operator list
   OpList["H"] = Hamiltonian;
   OpList["H_a"] = E;
   OpList["H_b"] = F;

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);

   pheap::Initialize(argv[5], 1, 65536, 655360);
   pheap::ShutdownPersistent(OList);
}
