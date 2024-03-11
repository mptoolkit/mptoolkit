// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// tensor/test/testadjoint.cpp
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

#include "tensor/tensor.h"
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

   // check that CH = adjoint(C)
   CHECK_CLOSE(adjoint(C), CH);

   // check that adjoint(CH) = -C
   CHECK_CLOSE(adjoint(CH), -C);

   // check that CH = -inv_adjoint(C)
   CHECK_CLOSE(inv_adjoint(C), -CH);

   // check that inv_adjoint(CH) = C
   CHECK_CLOSE(inv_adjoint(CH), C)(CH)(adjoint(CH))(inv_adjoint(CH))(C);


   // check that C * adoint(C) is what we expect (ie. N-2)
   IrredTensor<double> CCH(HubbardBasis, QN(0, 0));
   CCH(0,0)  = -2.0;
   CCH(1,1) = -1.0;

   CHECK_CLOSE(CCH, sqrt(2.0) * prod(C, adjoint(C), QN(0,0)));

   // check that adjoint(C) * C is what we expect (ie. N)
   IrredTensor<double> N(HubbardBasis, QN(0, 0));
   N(1,1)  = 1.0;
   N(2,2)  = 2.0;

   CHECK_CLOSE(N, sqrt(2.0) * prod(adjoint(C), C, QN(0,0)));

   // the 'inverse adjoint' should give the opposite sign for these spin-1/2 operators.
   CHECK_CLOSE(inv_adjoint(C), -CH);
   CHECK_CLOSE(CCH, -sqrt(2.0) * prod(C, inv_adjoint(C), QN(0,0)));
   CHECK_CLOSE(N, -sqrt(2.0) * prod(inv_adjoint(C), C, QN(0,0)));
}
