// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// tensor/test/testtensorprod.cpp
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

#include "tensor/tensorproduct.h"
#include "quantumnumbers/u1.h"
#include "quantumnumbers/su2.h"

using namespace Tensor;
using namespace QuantumNumbers;

int main()
{
   SymmetryList Symmetry("N:U(1),S:SU(2)");
   // helper function object for creating quantum numbers at compile time
   QNConstructor<U1,SU2> QN(Symmetry);

   BasisList HubbardBasis(Symmetry);
   HubbardBasis.push_back(QN(0, 0));   // empty
   HubbardBasis.push_back(QN(1, 0.5)); // single
   HubbardBasis.push_back(QN(2, 0));   // double

   // Define the annihilation operator
   IrredTensor<double> C(HubbardBasis, QN(-1, 0.5));
   C(0,1)    = sqrt(2.0);
   C(1,2)    = 1.0;

   // and the creation operator
   IrredTensor<double> CH(HubbardBasis, QN(1, 0.5));
   CH(1,0)    = 1.0;
   CH(2,1)    = -sqrt(2.0);

   IrredTensor<double> HI(HubbardBasis, QN(0, 0));
   HI(0,0) = 1.0;
   HI(1,1) = 1.0;
   HI(2,2) = 1.0;

   IrredTensor<double> HS(HubbardBasis, QN(0, 1));
   HS(1,1) = std::sqrt(0.75);

   BasisList SpinBasis(Symmetry);
   SpinBasis.push_back(QN(1, 0.5));

   ProductBasis<BasisList, BasisList> HSBasis(HubbardBasis, SpinBasis);
   CHECK_EQUAL(HSBasis.size(), 4);
   CHECK_EQUAL(HSBasis[0], QN(1, 0.5));
   CHECK_EQUAL(HSBasis[1], QN(2, 0));
   CHECK_EQUAL(HSBasis[2], QN(2, 1));
   CHECK_EQUAL(HSBasis[3], QN(3, 0.5));

   IrredTensor<double> SS(SpinBasis, QN(0, 1));
   SS(0,0) = sqrt(0.75);

   IrredTensor<double> SI(SpinBasis, QN(0, 0));
   SI(0,0) = 1.0;

   IrredTensor<double> HISS = tensor_prod(HI, SS);
   CHECK_EQUAL(HISS.Basis1(), HSBasis.Basis());
   CHECK_EQUAL(HISS.Basis2(), HSBasis.Basis());

   IrredTensor<double> HSSI = tensor_prod(HS, SI);
   CHECK_EQUAL(HSSI.Basis1(), HSBasis.Basis());
   CHECK_EQUAL(HSSI.Basis2(), HSBasis.Basis());

   IrredTensor<double> TotalS = HISS+HSSI;
   CHECK_CLOSE(TotalS(0,0), sqrt(0.75));
   CHECK_CLOSE(TotalS(2,2), sqrt(2.0));
   CHECK_CLOSE(TotalS(3,3), sqrt(0.75));
   CHECK_CLOSE(norm_frob_sq(TotalS), 9);  // 9 = 2 * (0.5*1.5) + 3 * (1*2) + 2 * (0.5*1.5)
}
