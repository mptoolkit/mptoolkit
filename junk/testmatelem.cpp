// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// junk/testmatelem.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "models/spin-su2.h"
#include "tensor/tensor_eigen.h"
#include "mps/state_component.h"
#include "tensor/regularize.h"

typedef IrredTensor<LinearAlgebra::Matrix<double>, VectorBasis, BasisList>
MixedOperator;

MatrixOperator Regularize(SimpleOperator const& Op)
{
   MixedOperator U = Tensor::Regularize(Op.Basis1());
   MixedOperator V = Tensor::Regularize(Op.Basis2());

   return triple_prod(U, Op, herm(V), QuantumNumber(Op.GetSymmetryList()), Op.TransformsAs());
}

int main()
{
   for (int s = 1; s <= 6; ++s)
      {

         double d = 2*s+1;

   SiteBlock B = CreateSU2SpinSite(s);

   QuantumNumbers::QNConstructor<QuantumNumbers::SU2> QN(B["I"].GetSymmetryList());

   SiteOperator Op = sqrt(3.0) * tensor_prod(B["S"], B["S"], QN(0));
   std::cout << Op << '\n';

   std::cout << Tensor::EigenvaluesHermitian(Regularize(Op)) << '\n';

}

}
