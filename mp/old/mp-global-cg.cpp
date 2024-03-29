// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/old/mp-global-cg.cpp
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

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp-algorithms/conjugategradient.h"
#include "mp/copyright.h"
#include <iostream>

struct MultTrunc
{
   typedef MPWavefunction const& argument_type;
   typedef MPWavefunction result_type;

   MultTrunc(MPOperator const& A, int NumStates) : A_(A), NumStates_(NumStates) {}

   result_type operator()(argument_type x) const
   {
      return prod(A_, x, x.TransformsAs(), NumStates_);
   }

   MPOperator const& A_;
   int NumStates_;
};

struct Inner
{
   std::complex<double> operator()(MPWavefunction const& x, MPWavefunction const& y) const
   {
      return overlap(x,y);
   }
};

struct Ident
{
   MPWavefunction const& operator()(MPWavefunction const& x) const
   {
      return x;
   }
};

int main(int argc, char** argv)
{
   if (argc < 5 || argc > 7)
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "usage: mp-global-cg <lattice> <operator> <psi> <rhs> [num-iter] [<maxstates>]\n";
      return 1;
   }

   TRACE(argc);

   pvalue_ptr<MPWavefunction> Psi = pheap::OpenPersistent(argv[3], 655360);
   pvalue_ptr<OperatorList> System = pheap::ImportHeap(argv[1]);
   std::string Operator = argv[2];
   pvalue_ptr<MPWavefunction> Rhs = pheap::ImportHeap(argv[4]);
   int NumIter = 4;
   if (argc >= 6) NumIter = boost::lexical_cast<int>(argv[5]);
   int MaxStates = DefaultMaxStates;
   if (argc == 7) MaxStates = boost::lexical_cast<int>(argv[6]);

   double tol = 1E-10;

   ConjugateGradient(*Psi.mutate(), MultTrunc((*System)[Operator], MaxStates), *Rhs, NumIter, tol,
                     Ident(), Inner(), Inner());

   std::cout.precision(14);
   std::cout << NumIter << " " << tol << '\n';

   pheap::ShutdownPersistent(Psi);
}
