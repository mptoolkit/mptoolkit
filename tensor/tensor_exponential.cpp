// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// tensor/tensor_exponential.cpp
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

#include "tensor_exponential.h"
#include "linearalgebra/exponential.h"
#include "tensor/regularize.h"

namespace Tensor
{

IrredTensor<std::complex<double>, BasisList, BasisList>
Exponentiate(IrredTensor<std::complex<double>, BasisList, BasisList> const& m)
{
   typedef IrredTensor<std::complex<double>, BasisList, BasisList> TensorType;
   PRECONDITION(is_scalar(m.TransformsAs()))("Can only exponentiate a scalar operator")(m.TransformsAs());
   PRECONDITION_EQUAL(m.Basis1(), m.Basis2());

   using QuantumNumbers::QuantumNumber;

   TensorType Result(m.Basis1(), m.Basis2());

   // enumerate the quantum numbers in m
   std::set<QuantumNumber> UsedQ = QuantumNumbersInBasis(m.Basis1());

   // linearize the basis
   for (std::set<QuantumNumber>::const_iterator Q = UsedQ.begin(); Q != UsedQ.end(); ++Q)
   {
      std::map<int, int> LinearMap = LinearizeQuantumNumberSubspace(m.Basis1(), *Q);
      LinearAlgebra::Matrix<std::complex<double> > M(LinearMap.size(), LinearMap.size(), 0.0);
      for (const_iterator<TensorType>::type I = iterate(m); I; ++I)
      {
         if (m.Basis1()[I.index()] != *Q)
            continue;
         for (const_inner_iterator<TensorType>::type J = iterate(I); J; ++J)
         {
            if (m.Basis2()[J.index2()] != *Q)
               continue;

            M(LinearMap[J.index1()], LinearMap[J.index2()]) = *J;
         }
      }

      M = LinearAlgebra::Exponentiate(1.0, M);
      for (std::map<int, int>::const_iterator I = LinearMap.begin(); I != LinearMap.end(); ++I)
      {
         for (std::map<int, int>::const_iterator J = LinearMap.begin(); J != LinearMap.end(); ++J)
         {
            std::complex<double> x = M(I->second, J->second);
            if (LinearAlgebra::norm_frob(x) > 1e-14)
               Result(I->first, J->first) = x;
         }
      }
   }

   return Result;
}

IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis>
exp(IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis> const& m)
{
   PRECONDITION(is_scalar(m.TransformsAs()))("Can only exponentiate a scalar operator")(m.TransformsAs());
   PRECONDITION_EQUAL(m.Basis1(), m.Basis2());

   if (!is_regular_basis(m.Basis1()))
   {
      IrredTensor<LinearAlgebra::Matrix<double>, VectorBasis, VectorBasis> X
         = Regularize(m.Basis1());

      return triple_prod(herm(X), exp(triple_prod(X, m, herm(X))), X);
   }

   IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis>
      Result = m;

   for (unsigned i = 0; i < Result.Basis1().size(); ++i)
   {
      if (!iterate_at(Result.data(), i,i))
      {
         int Dim = Result.Basis1().dim(i);
         // element is zero, so the exponential is the identity
         Result.data()(i,i) = LinearAlgebra::DiagonalMatrix<double>(Dim, Dim, 1.0);
      }
      else
      {
         Result.data()(i,i) = LinearAlgebra::Exponentiate(1.0, Result.data()(i,i));
      }
   }

   return Result;
}

} // namespace Tensor
