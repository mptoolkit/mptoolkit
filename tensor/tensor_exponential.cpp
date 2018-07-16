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
#include "blas/exponential.h"
#include "tensor/regularize.h"
#include "blas/functors.h"

namespace Tensor
{

IrredTensor<std::complex<double>, BasisList, BasisList>
exp(IrredTensor<std::complex<double>, BasisList, BasisList> const& m)
{
   typedef IrredTensor<std::complex<double>, BasisList, BasisList> TensorType;
   PRECONDITION(is_scalar(m.TransformsAs()))("Can only exponentiate a scalar operator")(m.TransformsAs());
   PRECONDITION_EQUAL(m.Basis1(), m.Basis2());

   using QuantumNumbers::QuantumNumber;

   TensorType Result(m.Basis1(), m.Basis2());

   // enumerate the quantum numbers in m
   std::set<QuantumNumber> UsedQ = QuantumNumbersInBasis(m.Basis1());

   // linearize the basis
   for (auto const& Q : UsedQ)
   {
      std::map<int, int> LinearMap = LinearizeQuantumNumberSubspace(m.Basis1(), Q);
      blas::Matrix<std::complex<double>> M(LinearMap.size(), LinearMap.size(), 0.0);
      for (auto const& I : m)
      {
         if (m.Basis1()[I.row()] != Q)
            continue;
	 for (auto const& J : I)
	 {
            if (m.Basis2()[J.col()] != Q)
               continue;

            M(LinearMap[I.row()], LinearMap[J.col()]) = J.value;
         }
      }

      M = exp(std::move(M));
      for (auto const& I : LinearMap)
      {
	 for (auto const& J : LinearMap)
	 {
            std::complex<double> x = M(I.second, J.second);
	    // TODO: culling matrix elements here needs some thought
            if (std::abs(x) > std::numeric_limits<double>::epsilon())
               Result.insert(I.first, J.first, x);
	 }
      }
   }

   return Result;
}

} // namespace Tensor
