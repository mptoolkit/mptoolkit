// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// tensor/test/testbasis.cpp
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

#include "tensor/basis.h"

#include "quantumnumbers/su2.h"
#include "quantumnumbers/u1.h"

using namespace QuantumNumbers;
using namespace Tensor;

int main()
{
   SymmetryList Q("S:SU(2)");
   QNConstructor<SU2> QN(Q);

   BasisList B(Q);
   B.push_back(QN(0.5));
   B.push_back(QN(0));
   B.push_back(QN(1));

   CHECK_EQUAL(B.size(), 3);
   CHECK_EQUAL(B[0], QN(0.5));
   CHECK_EQUAL(B[1], QN(0));
   CHECK_EQUAL(B[2], QN(1));


   std::list<QuantumNumber> ListTest;
   ListTest.push_back(QN(0.5));
   ListTest.push_back(QN(0));
   ListTest.push_back(QN(1));

   BasisList B2(ListTest.begin(), ListTest.end());

   CHECK_EQUAL(B, B2);

   BasisList Badj = adjoint(B);
   CHECK_EQUAL(B, Badj);

   BasisList Su2Vacuum = make_vacuum_basis(Q);
   CHECK_EQUAL(Su2Vacuum.size(), 1);
   CHECK_EQUAL(Su2Vacuum[0], QN(0));

   // tghe adjoint of U(1) quantum numbers is non-identity
   SymmetryList N("S:U(1)");
   QNConstructor<U1> NN(N);

   BasisList Bn(N);
   Bn.push_back(NN(-0.5));
   Bn.push_back(NN(0.5));
   Bn.push_back(NN(8));

   CHECK_EQUAL(Bn.size(), 3);
   CHECK_EQUAL(Bn[0], NN(-0.5));
   CHECK_EQUAL(Bn[1], NN(0.5));
   CHECK_EQUAL(Bn[2], NN(8));

   BasisList BnAdj = adjoint(Bn);
   CHECK_EQUAL(BnAdj.size(), 3);
   CHECK_EQUAL(BnAdj[0], NN(0.5));
   CHECK_EQUAL(BnAdj[1], NN(-0.5));
   CHECK_EQUAL(BnAdj[2], NN(-8));

   std::cout << B << show_projections(B);
   std::cout << Bn << show_projections(Bn);

   // VectorBasis

   std::vector<int> D(3);
   D[0] = 1;
   D[1] = 5;
   D[2] = 7;
   VectorBasis VB(B, D.begin(), D.end());

   CHECK_EQUAL(VB.size(), 3);
   CHECK_EQUAL(VB.total_dimension(), 13);
   CHECK_EQUAL(VB.total_degree(), 28);

   std::cout << VB << show_projections(VB);

   VectorBasis VBc(B);
   VBc.set_dim(0, 1);
   VBc.set_dim(1, 5);
   VBc.set_dim(2, 7);

   CHECK_EQUAL(VB, VBc);

}
