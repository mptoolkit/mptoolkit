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


template <typename T>
IrredTensor<T, VectorBasis, VectorBasis>
Regularize(VectorBasis const& b)
{
   typedef IrredTensor<T, VectorBasis, VectorBasis> ResultType;
   // Iterate through b and determine the total dimension of
   // each quantum number space, and also map the subspaces of b
   // onto a range of the total space.
   typedef std::map<QuantumNumbers::QuantumNumber, int> SizeMapType;
   SizeMapType SizeMap;
   std::vector<blas::Range> RangeOfSubspace;
   RangeOfSubspace.reserve(b.size());

   for (std::size_t i = 0; i < b.size(); ++i)
   {
      int Sz = SizeMap[b[i]];
      RangeOfSubspace.push_back(blas::range(Sz, Sz+b.dim(i)));
      SizeMap[b[i]] += b.dim(i);
   }

   // Now construct the basis
   SizeMapType IndexOfQ;
   {
      int i = 0;
      for (SizeMapType::const_iterator I = SizeMap.begin(); I != SizeMap.end(); ++I)
         IndexOfQ[I->first] = i++;
   }

   VectorBasis RegularBasis(SizeMap.begin(), SizeMap.end());
   ResultType Result(RegularBasis, b, QuantumNumbers::QuantumNumber(b.GetSymmetryList()));
   for (std::size_t i = 0; i < b.size(); ++i)
   {
      int Dest = IndexOfQ[b[i]];
      if (size1(Result(Dest, i)) == 0)
         Result(Dest, i) = T(RegularBasis.dim(Dest),
                             b.dim(i),
                             0.0);

      Result(Dest,i)(RangeOfSubspace[i], blas::all)
         = blas::identity_matrix<double>(b.dim(i));
   }

   return Result;
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

} // namespace Tensor
