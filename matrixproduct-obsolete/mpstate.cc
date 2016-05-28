// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// matrixproduct-obsolete/mpstate.cc
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

template <typename T>
BasicMPStateComponent<T>::BasicMPStateComponent(BasisList const& SBasis_, 
                                      VectorBasis const& V1, 
                                      VectorBasis const& V2)
   : SBasis(SBasis_), VBasis1(V1), VBasis2(V2), Data(SBasis_.size())
{
   for (std::size_t i = 0; i < SBasis.size(); ++i)
   {
      Data[i] = MatrixOperator(VBasis1, VBasis2, SBasis[i]);
   }
}

template <typename T>
BasicMPStateComponent<T> 
BasicMPStateComponent<T>::ConstructFullBasis1(BasisList const& S, VectorBasis const& Basis2)
{
   ProductBasis<BasisList, VectorBasis> FullLeftBasis(S, Basis2);
   QuantumNumbers::QuantumNumber Ident(S.GetSymmetryList());
   MPStateComponent Result(S, FullLeftBasis.Basis(), Basis2);
   for (std::size_t t = 0; t < FullLeftBasis.size(); ++t)
   {
      int s, b2;
      std::tie(s,b2) = FullLeftBasis.rmap(t);

      int Dim = FullLeftBasis.dim(t);
      CHECK_EQUAL(Dim, Basis2.dim(b2));

      // Make an identity matrix of the correct size
      LinearAlgebra::set_element(Result[s], t, b2, LinearAlgebra::identity_matrix<double>(Dim));
   }

   // check the normalization
   DEBUG_CHECK_CLOSE(norm_frob_sq(scalar_prod(Result, herm(Result))), 
                     FullLeftBasis.total_degree());

   return Result;
}

template <typename T>
BasicMPStateComponent<T> 
BasicMPStateComponent<T>::ConstructFullBasis2(VectorBasis const& Basis1, BasisList const& S)
{
   ProductBasis<VectorBasis, BasisList> FullRightBasis(Basis1, adjoint(S));
   MPStateComponent Result(S, Basis1, FullRightBasis.Basis());
   for (std::size_t t = 0; t < FullRightBasis.size(); ++t)
   {
      int s, b1;
      std::tie(b1,s) = FullRightBasis.rmap(t);

      int Dim = FullRightBasis.dim(t);
      CHECK_EQUAL(Dim, Basis1.dim(b1));

      // Make an identity matrix of the correct size
      set_element(Result[s], b1, t, 
                                 std::sqrt(double(degree(FullRightBasis[t])) / degree(Basis1[b1])) 
                                 * LinearAlgebra::identity_matrix<double>(Dim));
   }

   // check the normalization
   DEBUG_CHECK_CLOSE(norm_frob_sq(scalar_prod(herm(Result), Result)), 
                     FullRightBasis.total_degree());

   return Result;
}

template <typename T>
PStream::opstream& operator<<(PStream::opstream& out, BasicMPStateComponent<T> const& Op)
{
   return out << Op.SBasis << Op.VBasis1 << Op.VBasis2 << Op.Data;
}

template <typename T>
PStream::ipstream& operator>>(PStream::ipstream& in, BasicMPStateComponent<T>& Op)
{
   return in >> Op.SBasis >> Op.VBasis1 >> Op.VBasis2 >> Op.Data;
}
