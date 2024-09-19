// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// attic/matrixproduct-obsolete/test/testdensity.cpp
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

#include "matrixproduct/density.h"
#include "quantumnumbers/su2.h"

int main()
{
   SymmetryList Symmetry("S:SU(2)");
   QuantumNumbers::QNConstructor<QuantumNumbers::SU2> QN(Symmetry);
   BasisList B(Symmetry);
   B.push_back(QN(1));
   B.push_back(QN(1));

   SimpleOperator M(B, B, QN(0));
   M(0,0) = 0;
   M(0,1) = 0.5;
   M(1,0) = 0.5;
   M(1,1) = 0;

   DensityMatrix<SimpleOperator> DM(M);

   DM.DensityMatrixReport(std::cout);
}
