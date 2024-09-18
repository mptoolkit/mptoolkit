// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// attic/matrixproduct-obsolete/test/testscalardirectprod.cpp
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

#include "matrixproduct/mpstate.h"
#include "quantumnumbers/su2.h"

int main()
{
   // construct an AKLT state
   SymmetryList Symmetry("S:SU(2)");
   QuantumNumbers::QNConstructor<QuantumNumbers::SU2> QN(Symmetry);

   BasisList LocalBasis(Symmetry);
   LocalBasis.push_back(QN(1));  // S=1 chain

   VectorBasis MatrixBasis(Symmetry);
   MatrixBasis.push_back(QN(0.5));    // spin 1/2 edge states

   MPStateComponent AKLT(LocalBasis, MatrixBasis, MatrixBasis);
   AKLT[0] = MatrixOperator(MatrixBasis, QN(1));
   AKLT[0](0,0) = LinearAlgebra::Matrix<double>(1,1);
   AKLT[0](0,0)(0,0) = sqrt(2.0); // sqrt(s*(s+1)) = sqrt(2.0)

   TRACE(LinearAlgebra::scalar_direct_prod(herm(AKLT), AKLT));
   TRACE(LinearAlgebra::scalar_direct_prod(AKLT, herm(AKLT)));
}
