// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mps/state_component.cc
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

template <typename T>
BasicStateComponent<T>::BasicStateComponent(BasisList const& SBasis_,
                                      VectorBasis const& V1,
                                      VectorBasis const& V2)
   : SBasis(SBasis_), VBasis1(V1), VBasis2(V2)
{
   for (std::size_t i = 0; i < SBasis.size(); ++i)
   {
      Data.push_back(value_type(VBasis1, VBasis2, SBasis[i]));
   }
}

template <typename T>
BasicStateComponent<T>::BasicStateComponent(BasisList const& SBasis_,
                                      VectorBasis const& V1,
					    VectorBasis const& V2, DataType const& D)
   : SBasis(SBasis_), VBasis1(V1), VBasis2(V2)
{
   for (std::size_t i = 0; i < SBasis.size(); ++i)
   {
      Data.push_back(copy(D[i]));
   }
}

template <typename T>
BasicStateComponent<T>
BasicStateComponent<T>::ConstructFullBasis1(BasisList const& S, VectorBasis const& Basis2)
{
   ProductBasis<BasisList, VectorBasis> FullLeftBasis(S, Basis2);
   QuantumNumbers::QuantumNumber Ident(S.GetSymmetryList());
   StateComponent Result(S, FullLeftBasis.Basis(), Basis2);
   for (std::size_t t = 0; t < FullLeftBasis.size(); ++t)
   {
      int s, b2;
      std::tie(s,b2) = FullLeftBasis.rmap(t);

      int Dim = FullLeftBasis.dim(t);
      CHECK_EQUAL(Dim, Basis2.dim(b2));

      Result[s].insert(t, b2, value_type::make_identity(Dim));
   }

   // check the normalization
   //   DEBUG_CHECK_CLOSE(norm_frob_sq(scalar_prod(Result, herm(Result))),
   //                     FullLeftBasis.total_degree());
   //DEBUG_CHECK_CLOSE(norm_frob_sq(Result), FullLeftBasis.total_degree());

   return Result;
}

template <typename T>
void
BasicStateComponent<T>::delta_shift(QuantumNumber const& q)
{
   CHECK_EQUAL(degree(q), 1);
   VBasis1.delta_shift(q);
   VBasis2.delta_shift(q);
   for (auto& c : Data)
   {
      c.delta_shift(q);
   }
}

template <typename T>
BasicStateComponent<T>
BasicStateComponent<T>::ConstructFullBasis2(VectorBasis const& Basis1, BasisList const& S)
{
   ProductBasis<VectorBasis, BasisList> FullRightBasis(Basis1, adjoint(S));
   StateComponent Result(S, Basis1, FullRightBasis.Basis());
   for (std::size_t t = 0; t < FullRightBasis.size(); ++t)
   {
      int s, b1;
      std::tie(b1,s) = FullRightBasis.rmap(t);

      int Dim = FullRightBasis.dim(t);
      CHECK_EQUAL(Dim, Basis1.dim(b1));

      // Make an identity matrix of the correct size
      Result[s].insert(b1, t,
                       std::sqrt(double(degree(FullRightBasis[t])) / degree(Basis1[b1]))
                       * value_type::make_identity(Dim));
   }

   // check the normalization
   //DEBUG_CHECK_CLOSE(norm_frob_sq(scalar_prod(herm(Result), Result)),
   //                FullRightBasis.total_degree());

   return Result;
}

template <typename T>
PStream::opstream& operator<<(PStream::opstream& out, BasicStateComponent<T> const& Op)
{
   return out << Op.SBasis << Op.VBasis1 << Op.VBasis2 << Op.Data;
}

template <typename T>
PStream::ipstream& operator>>(PStream::ipstream& in, BasicStateComponent<T>& Op)
{
   return in >> Op.SBasis >> Op.VBasis1 >> Op.VBasis2 >> Op.Data;
}

#if 0
namespace LinearAlgebra
{

template <typename T, typename U, typename Func>
MatrixOperator
ScalarProd<BasicStateComponent<T>, HermitianProxy<BasicStateComponent<U>>, Func>::
operator()(BasicStateComponent<T> const& A, HermitianProxy<BasicStateComponent<U>> const& B) const
{
   DEBUG_PRECONDITION_EQUAL(A.LocalBasis(), B.base().LocalBasis());
   DEBUG_PRECONDITION_EQUAL(A.Basis2(), B.base().Basis2());

   QuantumNumbers::QuantumNumber Ident(A.GetSymmetryList());
   MatrixOperator Result(A.Basis1(), B.base().Basis1(), Ident);

   for (std::size_t s = 0; s < A.LocalBasis().size(); ++s)
   {
      Result += scalar_prod(A[s], herm(B.base()[s]));
   }
   return Result;
}



} //namespace LinearAlgebra
#endif
