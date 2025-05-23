// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// tensor/tensor_eigen.cc
//
// Copyright (C) 2004-2024 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2022 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#include "linearalgebra/matrix_utility.h"
#include "regularize.h"

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
      if (A(i,i).size() > 0)
         Result(i,i).diagonal() = real_diagonal_vector(A(i,i));
      else
         Result(i,i).diagonal() = LinearAlgebra::Vector<T>(A.Basis1().dim(i), T{});
   }
   Result.debug_check_structure();
   return Result;
}

template <typename T>
void
OrthogonalizeRowsAgainstRegular(IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis>& X, IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis> const& Y)
{
   DEBUG_CHECK(is_regular_basis(X.Basis1()));
   DEBUG_CHECK(is_regular_basis(X.Basis2()));
   DEBUG_CHECK(is_regular_basis(Y.Basis1()));
   DEBUG_CHECK_EQUAL(X.Basis2(), Y.Basis2());
   for (auto I = iterate(X); I; ++I)
   {
      for (auto J = iterate(I); J; ++J)
      {
         // see if any rows in Y have an entry at J.index2()
         bool Found = false;
         for (int i = 0; i < Y.Basis1().size(); ++i)
         {
            auto K = iterate_at(Y.data(), i, J.index2());
            if (K)
            {
               OrthogonalizeRowsAgainst(*J, *K);
               Found = true;
            }
         }
         if (!Found)
         {
            *J = LQ_FactorizeThin(std::move(*J)).second;
         }
      }
   }
}

template <typename T>
void
OrthogonalizeRowsAgainst(IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis>& X, IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis> const& Y)
{
   CHECK(is_scalar(X.TransformsAs()));
   CHECK(is_scalar(Y.TransformsAs()));
   DEBUG_CHECK_EQUAL(X.Basis2(), Y.Basis2());
   if (is_regular_basis(X.Basis1()) && is_regular_basis(X.Basis2()) && is_regular_basis(Y.Basis1()))
   {
      OrthogonalizeRowsAgainstRegular(X, Y);
      return;
   }
   // else

   Regularizer R1(X.Basis1());
   Regularizer R2(X.Basis2());
   Regularizer S1(Y.Basis1());
   auto XX = RegularizeBasis12(R1, std::move(X), R2);
   auto YY = RegularizeBasis12(S1, Y, R2);
   OrthogonalizeRowsAgainstRegular(XX, YY);
   X = UnregularizeBasis12(R1, XX, R2);
}

template <typename T>
void
OrthogonalizeColsAgainstRegular(IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis>& X, IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis> const& Y)
{
   DEBUG_CHECK(is_regular_basis(X.Basis1()));
   DEBUG_CHECK(is_regular_basis(X.Basis2()));
   DEBUG_CHECK(is_regular_basis(Y.Basis2()));
   DEBUG_CHECK_EQUAL(X.Basis1(), Y.Basis1());
   for (auto I = iterate(X); I; ++I)
   {
      for (auto J = iterate(I); J; ++J)
      {
         auto K = Y.data().vec()[J.index1()].iterate();
         if (K)
         {
            OrthogonalizeColsAgainst(*J, *K);
         }
         else
         {
            *J = QR_FactorizeThin(std::move(*J)).first;
         }
      }
   }
}

template <typename T>
void
OrthogonalizeColsAgainst(IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis>& X, IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis> const& Y)
{
   CHECK(is_scalar(X.TransformsAs()));
   CHECK(is_scalar(Y.TransformsAs()));
   DEBUG_CHECK_EQUAL(X.Basis1(), Y.Basis1());
   if (is_regular_basis(X.Basis1()) && is_regular_basis(X.Basis2()) && is_regular_basis(Y.Basis2()))
   {
      OrthogonalizeColsAgainstRegular(X, Y);
      return;
   }
   // else

   Regularizer R1(X.Basis1());
   Regularizer R2(X.Basis2());
   Regularizer S2(Y.Basis2());
   auto XX = RegularizeBasis12(R1, std::move(X), R2);
   auto YY = RegularizeBasis12(R1, Y, S2);
   OrthogonalizeColsAgainstRegular(XX, YY);
   X = UnregularizeBasis12(R1, XX, R2);
}

} // namespace Tensor
