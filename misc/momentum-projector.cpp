// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// misc/momentum-projector.cpp
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
#include "matrixproduct/operatoratsite.h"
#include "pheap/pheap.h"

typedef std::complex<double> complex;
typedef LinearAlgebra::Matrix<complex> ComplexMatrix;

ComplexMatrix Hamiltonian(double t, std::vector<double> const& U, bool Periodic = false)
{
   int size = U.size();
   ComplexMatrix Result(size, size, 0.0);
   for (int i = 0; i < size-1; ++i)
   {
      Result(i,i+1) = Result(i+1,i) = -t;
   }
   for (int i = 0; i < size; ++i)
   {
      Result(i,i) = U[i];
   }
   if (Periodic)
      Result(0,size-1) = Result(size-1,0) = -t;

   return Result;
}

int main(int argc, char** argv)
{
   if (argc != 2)
   {
      std::cerr << "usage: projector <lattice>\n";
      return 1;
   }

   // load the lattice file
   pvalue_ptr<OperatorList> pSystem = pheap::OpenPersistent(argv[1], 655360);
   OperatorList System = *pSystem;

   // get the lattice size
   int LatticeSize = System.size();

   std::vector<double> U(LatticeSize, 0);
   ComplexMatrix H = Hamiltonian(1.0, U);

   OperatorAtSite<OperatorList const> CH(System, "CH");
   OperatorAtSite<OperatorList const> C(System, "C");

   ComplexMatrix Diag = H;
   LinearAlgebra::Vector<double> EVal = DiagonalizeHermitian(Diag);
   TRACE(EVal);

   QuantumNumber Ident(System.GetSymmetryList());

   // Add the momentum space operators
   OperatorAtSite<OperatorList> CkH(System, "CkH");
   OperatorAtSite<OperatorList> Ck(System, "Ck");
   OperatorAtSite<OperatorList> Nk(System, "Nk");
   for (int i = 0; i < LatticeSize; ++i)
   {
      std::cout << "Working..." << std::endl;
      Ck(i) = MPOperator();
      CkH(i) = MPOperator();
      LinearAlgebra::Vector<std::complex<double> > v = Diag(i,LinearAlgebra::all);
      for (int j = 0; j < LatticeSize; ++j)
      {
         Ck(i) += v[j] * C(j+1);
         CkH(i) += conj(v[j]) * CH(j+1);
      }
      Nk(i) = prod(CkH(i), Ck(i), Ident);
   }

   // save the lattice
   *pSystem.mutate() = System;
   pheap::ShutdownPersistent(pSystem);
}
