// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// tensor/regularize.cc
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

namespace Tensor
{

inline
int Regularizer::IndexOf(int i) const
{
   return BasisMappingIndex[i];
}

inline
LinearAlgebra::Range Regularizer::RangeOf(int i) const
{
   return BasisMappingRange[i];
}

template <typename T>
IrredTensor<T, BasisList, BasisList>
map_1x1_operator(IrredTensor<LinearAlgebra::Matrix<T>, BasisList, BasisList> const& Op)
{
   typedef IrredTensor<LinearAlgebra::Matrix<T>, BasisList, BasisList> MatOpType;
   IrredTensor<T, BasisList, BasisList> Result(Op.Basis1(), Op.Basis2(), Op.TransformsAs());

   for (typename const_iterator<MatOpType>::type I = iterate(Op); I; ++I)
   {
      for (typename const_inner_iterator<MatOpType>::type J = iterate(I); J; ++J)
      {
         DEBUG_PRECONDITION_EQUAL(size1(*J), 1);
         DEBUG_PRECONDITION_EQUAL(size2(*J), 1);
         Result(J.index1(), J.index2()) = (*J)(0,0);
      }
   }
   return Result;
}

template <typename T>
IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis>
ToMatrixOperator(IrredTensor<T, BasisList, BasisList> const& Op)
{
   IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis>
      Result(VectorBasis(Op.Basis1()), VectorBasis(Op.Basis2()), Op.TransformsAs());
   Result.data() = LinearAlgebra::scalar(LinearAlgebra::Matrix<T>(1,1,1.0)) * Op.data();
   return Result;
}

template <typename T>
IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis>
RegularizeBasis1(Regularizer const& R, IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis> const& M)
{
   IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis> Result(R.Basis(), M.Basis2(), M.TransformsAs());
   for (auto I = iterate(M); I; ++I)
   {
      for (auto J = iterate(I); J; ++J)
      {
         auto r = iterate_at(Result.data(), R.IndexOf(J.index1()), J.index2());
         if (!r)
         {
            Result(R.IndexOf(J.index1()), J.index2()) = LinearAlgebra::Matrix<T>(R.Basis().dim(R.IndexOf(J.index1())), M.Basis2().dim(J.index2()), 0.0);
            r = iterate_at(Result.data(), R.IndexOf(J.index1()), J.index2());
         }
         (*r)(R.RangeOf(J.index1()), LinearAlgebra::all) = *J;
      }
   }
   return Result;
}

template <typename T>
IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis>
UnregularizeBasis1(Regularizer const& R, IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis> const& M)
{
   IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis> Result(R.OriginalBasis(), M.Basis2(), M.TransformsAs());
   for (int i = 0; i < Result.Basis1().size(); ++i)
   {
      for (int j = 0; j < Result.Basis2().size(); ++j)
      {
         auto r = iterate_at(M.data(), R.IndexOf(i), j);
         if (r)
         {
            Result(i,j) = (*r)(R.RangeOf(i), LinearAlgebra::all);
         }
      }
   }
   return Result;
}

template <typename T>
IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis>
RegularizeBasis2(IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis> const& M, Regularizer const& R)
{
   IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis> Result(M.Basis1(), R.Basis(), M.TransformsAs());
   for (auto I = iterate(M); I; ++I)
   {
      for (auto J = iterate(I); J; ++J)
      {
         auto r = iterate_at(Result.data(), J.index1(), R.IndexOf(J.index2()));
         if (!r)
         {
            Result(J.index1(), R.IndexOf(J.index2())) = LinearAlgebra::Matrix<T>(M.Basis1().dim(J.index1()), R.Basis().dim(R.IndexOf(J.index2())), 0.0);
            r = iterate_at(Result.data(), J.index1(), R.IndexOf(J.index2()));
         }
         (*r)(LinearAlgebra::all, R.RangeOf(J.index2())) = *J;
      }
   }
   return Result;
}

template <typename T>
IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis>
UnregularizeBasis2(IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis> const& M, Regularizer const& R)
{
   IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis> Result(M.Basis1(), R.OriginalBasis(), M.TransformsAs());
   for (int i = 0; i < Result.Basis1().size(); ++i)
   {
      for (int j = 0; j < Result.Basis2().size(); ++j)
      {
         auto r = iterate_at(M.data(), i, R.IndexOf(j));
         if (r)
         {
            Result(i,j) = (*r)(LinearAlgebra::all, R.RangeOf(j));
         }
      }
   }
   return Result;
}

template <typename T>
IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis>
RegularizeBasis12(Regularizer const& R1, IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis> const& M, Regularizer const& R2)
{
   IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis> Result(R1.Basis(), R2.Basis(), M.TransformsAs());
   for (auto I = iterate(M); I; ++I)
   {
      for (auto J = iterate(I); J; ++J)
      {
         auto r = iterate_at(Result.data(), R1.IndexOf(J.index1()), R2.IndexOf(J.index2()));
         if (!r)
         {
            Result(R1.IndexOf(J.index1()), R2.IndexOf(J.index2())) = LinearAlgebra::Matrix<T>(R1.Basis().dim(R1.IndexOf(J.index1())), R2.Basis().dim(R2.IndexOf(J.index2())), 0.0);
            r = iterate_at(Result.data(), R1.IndexOf(J.index1()), R2.IndexOf(J.index2()));
         }
         (*r)(R1.RangeOf(J.index1()), R2.RangeOf(J.index2())) = *J;
      }
   }
   return Result;
}

template <typename T>
IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis>
UnregularizeBasis12(Regularizer const& R1, IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis> const& M, Regularizer const& R2)
{
   IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis> Result(R1.OriginalBasis(), R2.OriginalBasis(), M.TransformsAs());
   for (int i = 0; i < Result.Basis1().size(); ++i)
   {
      for (int j = 0; j < Result.Basis2().size(); ++j)
      {
         auto r = iterate_at(M.data(), R1.IndexOf(i), R2.IndexOf(i));
         if (r)
         {
            Result(i,j) = (*r)(R1.RangeOf(i), R2.RangeOf(j));
         }
      }
   }
   return Result;
}


template <typename T>
IrredTensor<LinearAlgebra::DiagonalMatrix<T>, VectorBasis, VectorBasis, Tensor::DiagonalStructure>
RegularizeBasis12(Regularizer const& R1, IrredTensor<LinearAlgebra::DiagonalMatrix<T>, VectorBasis, VectorBasis, Tensor::DiagonalStructure> const& M, Regularizer const& R2)
{
   PANIC("Not implemented");
}

template <typename T>
IrredTensor<LinearAlgebra::DiagonalMatrix<T>, VectorBasis, VectorBasis, Tensor::DiagonalStructure>
UnregularizeBasis12(Regularizer const& R1, IrredTensor<LinearAlgebra::DiagonalMatrix<T>, VectorBasis, VectorBasis, Tensor::DiagonalStructure> const& M, Regularizer const& R2)
{
   PANIC("Not implemented");
}

} // namespace Tensor
