// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// matrixproduct-obsolete/test/testpackunpack.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "matrixproduct/packunpack.h"
#include "quantumnumbers/su2.h"
#include "matrixproduct/wavefunc-utils.h"

int main()
{
   SymmetryList Symmetry("S:SU(2)");
   QuantumNumbers::QNConstructor<QuantumNumbers::SU2> QN(Symmetry);
   VectorBasis B1(Symmetry);
   B1.push_back(QN(0), 2);
   B1.push_back(QN(0.5), 4);
   B1.push_back(QN(1), 5);
   B1.push_back(QN(2), 7);
   B1.push_back(QN(2.5), 8);

   VectorBasis B2(Symmetry);
   B2.push_back(QN(0), 3);
   B2.push_back(QN(1.5), 8);
   B2.push_back(QN(1), 2);
   B2.push_back(QN(2), 4);
   B2.push_back(QN(2.5), 7);

   MatrixOperator M = MakeRandomMatrixOperator(B1, B2, QN(1));

   PackMatrixOperator Packer(M);

   LinearAlgebra::Vector<std::complex<double> > v(Packer.size());
   Packer.pack(M, data(v));

   CHECK_CLOSE(norm_frob(M), norm_frob(v));

   MatrixOperator MCheck = Packer.unpack(data(v));
   CHECK_CLOSE(norm_frob(M-MCheck), 1e-14);
}
