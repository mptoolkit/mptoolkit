// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// tensor/tensor_eigen.cc
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

#include "linearalgebra/matrix_utility.h"

namespace Tensor
{

template <typename T>
IrredTensor<LinearAlgebra::DiagonalMatrix<T>, VectorBasis, VectorBasis>
InvertDiagonal(IrredTensor<LinearAlgebra::DiagonalMatrix<T>,
               VectorBasis, VectorBasis> const& Op, double Cutoff)
{
   IrredTensor<LinearAlgebra::DiagonalMatrix<T>, VectorBasis, VectorBasis>
      Result(Op.Basis2(), Op.Basis1());

   for (unsigned i = 0; i < Op.Basis1().size(); ++i)
   {
      if (iterate_at(Op.data(), i, i))
      {
         Result(i,i) = LinearAlgebra::DiagonalMatrix<T>(Op.Basis2().dim(i), Op.Basis1().dim(i), 0.0);
         for (int j = 0; j < std::min(Op.Basis1().dim(i), Op.Basis2().dim(i)); ++j)
         {
            T x = Op(i,i)(j,j);
            if (LinearAlgebra::norm_frob(x) < Cutoff)
               x = 0;
            else
               x = 1.0 / x;
            Result(i,i)(j,j) = x;
         }
      }
   }
   return Result;
}

template <typename T>
IrredTensor<LinearAlgebra::DiagonalMatrix<T>, VectorBasis, VectorBasis>
InvertDiagonal(IrredTensor<LinearAlgebra::DiagonalMatrix<T>,
               VectorBasis, VectorBasis> const& Op)
{
   return InvertDiagonal(Op, std::sqrt(std::numeric_limits<double>::epsilon()));
}

template <typename T>
IrredTensor<LinearAlgebra::DiagonalMatrix<T>, VectorBasis, VectorBasis, Tensor::DiagonalStructure>
ExtractDiagonal(IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis> A)
{
   DEBUG_CHECK_EQUAL(A.Basis1(), A.Basis2());
   CHECK(is_scalar(A.TransformsAs()));
   IrredTensor<LinearAlgebra::DiagonalMatrix<T>, VectorBasis, VectorBasis, Tensor::DiagonalStructure> Result(A.Basis1(), A.Basis2());
   for (int i = 0; i < A.Basis1().size(); ++i)
   {
      Result(i,i).diagonal() = diagonal_vector(A(i,i));
   }
   return Result;
}

template <typename T>
IrredTensor<LinearAlgebra::DiagonalMatrix<T>, VectorBasis, VectorBasis, Tensor::DiagonalStructure>
ExtractRealDiagonal(IrredTensor<LinearAlgebra::Matrix<std::complex<T>>, VectorBasis, VectorBasis> A)
{
   DEBUG_CHECK_EQUAL(A.Basis1(), A.Basis2());
   CHECK(is_scalar(A.TransformsAs()));
   IrredTensor<LinearAlgebra::DiagonalMatrix<T>, VectorBasis, VectorBasis, Tensor::DiagonalStructure> Result(A.Basis1(), A.Basis2());
   for (int i = 0; i < A.Basis1().size(); ++i)
   {
      Result(i,i).diagonal() = real_diagonal_vector(A(i,i));
   }
   return Result;
}

} // namespace Tensor
