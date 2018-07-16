// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mps/packunpack.cpp
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

#include "packunpack.h"
#include <cstring>

//
// PackMatrixOperator
//

void PackMatrixOperator::Initialize(VectorBasis const& Basis1,
                                     VectorBasis const& Basis2,
                                     QuantumNumbers::QuantumNumber const& q)
{
   B1_ = Basis1;
   B2_ = Basis2;
   q_ = q;
   OffsetMatrix_ = blas::Matrix<int>(Basis1.size(), Basis2.size(), -1);
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
   this->Initialize(Basis1, Basis2, q);
}

PackMatrixOperator::PackMatrixOperator(MatrixOperator const& m)
{
   Initialize(m.Basis1(), m.Basis2(), m.TransformsAs());
}

void
PackMatrixOperator::pack(cpu::Matrix const& m, value_type* Iter, double Factor)
{
   // pack in column-major format
   for (int i = 0; i < m.cols(); ++i)
   {
      vector_view(m.rows(), 1, Iter+i*m.leading_dimension())
	 = m.column(i);
   }
}

void
PackMatrixOperator::unpack(cpu::Matrix& m, value_type const* Iter, double Factor)
{
   // pack in column-major format
   for (int i = 0; i < m.cols(); ++i)
   {
      vector_view(m.rows(), 1, Iter+i*m.leading_dimension())
	 = Factor * m.column(i);
   }
}

void
PackMatrixOperator::unpack(cpu::Matrix& m, value_type const* Iter, double Factor)
{
   // pack in column-major format
   for (int i = 0; i < m.cols(); ++i)
   {
      m.column(i) = Factor * const_vector_view(m.rows(), 1, Iter+i*m.leading_dimension());
   }
}

PackMatrixOperator::value_type*
PackMatrixOperator::pack(MatrixOperator const& m, value_type* Iter) const
{
   for (auto const& I : OffsetArray_)
   {
      auto J = m.row(I.r).find(I.c);
      if (J != m.row(I.r).end())
      {
      {
	 DEBUG_CHECK_EQUAL(I.Size, J.value.rows()*J.value.cols());
         // normalization correction is sqrt(degree(B1[I->r]))
	 pack(get_wait(J.value), Iter+I.offset, std::sqrt(double(degree(B1_[I.r]))));
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
   MatrixOperator Result(B1_, B2_, q_);
   for (auto const& I : OffsetArray_)
   {
      blas::Matrix<value_type> m(B1_.dim(I->r), B2_.dim(I->c));
      unpack(m, Iter+I.Offset, 1.0 / std::sqrt(double(degree(B1_[I->r]))));
      Result.set(I.r, I.c, std::move(m));
   }
   return Result;
}

//
// PackStateComponent
//

void PackStateComponent::Initialize(BasisList const& LocalBasis, VectorBasis const& Basis1,
                                    VectorBasis const& Basis2)
{
   B_ = LocalBasis;
   B1_ = Basis1;
   B2_ = Basis2;
   OffsetMatrix_ = OffsetMatrixType(B_.size(),
                                    blas::Matrix<int>(Basis1.size(), Basis2.size(), -1));
   unsigned Offset = 0;
   for (unsigned q = 0; q < B_.size(); ++q)
   {
      for (unsigned i = 0; i < Basis1.size(); ++i)
      {
         for (unsigned j = 0; j < Basis2.size(); ++j)
         {
            if (is_transform_target(Basis2[j], B_[q], Basis1[i]))
            {
               unsigned Sz = Basis1.dim(i) * Basis2.dim(j);
               OffsetMatrix_[q](i,j) = Offset;
               OffsetArray_.push_back(OffsetRecType(q,i,j,Offset,Sz));
               Offset += Sz;
            }
         }
      }
   }
   Size_ = Offset;
}

PackStateComponent::PackStateComponent(BasisList const& LocalBasis, VectorBasis const& Basis1,
                                       VectorBasis const& Basis2)
{
   this->Initialize(LocalBasis, Basis1, Basis2);
}

PackStateComponent::PackStateComponent(StateComponent const& m)
{
   Initialize(m.LocalBasis(), m.Basis1(), m.Basis2());
}

PackStateComponent::value_type*
PackStateComponent::pack(StateComponent const& m, value_type* Iter) const
{
   DEBUG_CHECK_EQUAL(m.LocalBasis(), B_);
   DEBUG_CHECK_EQUAL(m.Basis1(), B1_);
   DEBUG_CHECK_EQUAL(m.Basis2(), B2_);
   typedef LinearAlgebra::VectorMemProxy<value_type> VecProxy;
   typedef LinearAlgebra::VectorMemProxy<value_type const> ConstVecProxy;
   for (auto const& I : OffsetArray_)
   {
      const_inner_iterator<MatrixOperator>::type J = iterate_at(m[I->q].data(), I->r, I->c);
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

StateComponent
PackStateComponent::unpack(value_type const* Iter) const
{
   //   typedef LinearAlgebra::VectorMemProxy<value_type> VecProxy;
   StateComponent Result(B_, B1_, B2_);
   for (OffsetArrayType::const_iterator I = OffsetArray_.begin();
        I != OffsetArray_.end(); ++I)
   {
      LinearAlgebra::Matrix<value_type> m(B1_.dim(I->r), B2_.dim(I->c));
      std::memcpy(data(m), Iter+I->Offset, I->Size * sizeof(value_type));
      // reverse the normalization correction
      Result[I->q](I->r, I->c) = (1.0 / std::sqrt(double(degree(B1_[I->r])))) * m;
   }
   return Result;
}
