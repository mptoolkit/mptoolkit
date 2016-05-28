// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// matrixproduct-obsolete/packunpack.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "packunpack.h"
#include <cstring>

void PackMatrixOperator:: Initialize(VectorBasis const& Basis1, 
                                     VectorBasis const& Basis2, 
                                     QuantumNumbers::QuantumNumber const& q)
{
   B1_ = Basis1;
   B2_ = Basis2;
   q_ = q;
   OffsetMatrix_ = LinearAlgebra::Matrix<int>(Basis1.size(), Basis2.size(), -1);
   unsigned Offset = 0;
   for (unsigned i = 0; i < Basis1.size(); ++i)
   {
      for (unsigned j = 0; j < Basis2.size(); ++j)
      {
         if (is_transform_target(Basis2[j], q, Basis1[i]))
         {
            unsigned Sz = Basis1.dim(i) * Basis2.dim(j);
            OffsetMatrix_(i,j) = Offset;
            OffsetArray_.push_back(OffsetRecType(i,j,Offset,Sz));
            Offset += Sz;
         }
      }
   }
   Size_ = Offset;
}

PackMatrixOperator::PackMatrixOperator(VectorBasis const& Basis1, 
                                       VectorBasis const& Basis2, 
                                       QuantumNumbers::QuantumNumber const& q)
{
   Initialize(Basis1, Basis2, q);
}

PackMatrixOperator::PackMatrixOperator(MatrixOperator const& m)
{
   Initialize(m.Basis1(), m.Basis2(), m.TransformsAs());
}

PackMatrixOperator::value_type* 
PackMatrixOperator::pack(MatrixOperator const& m, value_type* Iter) const
{
   typedef LinearAlgebra::VectorMemProxy<value_type> VecProxy;
   typedef LinearAlgebra::VectorMemProxy<value_type const> ConstVecProxy;
   for (OffsetArrayType::const_iterator I = OffsetArray_.begin();
        I != OffsetArray_.end(); ++I)
   {
      const_inner_iterator<MatrixOperator>::type J = iterate_at(m.data(), I->r, I->c);
      if (J)
      {
         DEBUG_CHECK_EQUAL(I->Size, size1(*J)*size2(*J));
         // normalization correction is sqrt(degree(B1[I->r]))
         VecProxy(Iter+I->Offset, I->Size) 
            = std::sqrt(double(degree(B1_[I->r]))) * ConstVecProxy(data(*J), I->Size);
      }
      else
      {
         std::memset(Iter+I->Offset, 0, I->Size * sizeof(value_type));
      }
   }
   return Iter+Size_;
}

MatrixOperator
PackMatrixOperator::unpack(value_type const* Iter) const
{
   //   typedef LinearAlgebra::VectorMemProxy<value_type> VecProxy;
   MatrixOperator Result(B1_, B2_, q_);
   for (OffsetArrayType::const_iterator I = OffsetArray_.begin();
        I != OffsetArray_.end(); ++I)
   {
      LinearAlgebra::Matrix<value_type> m(B1_.dim(I->r), B2_.dim(I->c));
      std::memcpy(data(m), Iter+I->Offset, I->Size * sizeof(value_type));
      // reverse the normalization correction
      Result(I->r, I->c) = (1.0 / std::sqrt(double(degree(B1_[I->r])))) * m;
   }
   return Result;
}
