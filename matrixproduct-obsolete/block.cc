// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// matrixproduct-obsolete/block.cc
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "siteoperator/siteoperator.h"

template <typename T>
void CoerceSymmetryList(Block<T>& block, QuantumNumbers::SymmetryList const& sl)
{
   block.CoerceSymmetryList(sl);
}

template <typename T>
void
Block<T>::CoerceSymmetryList(QuantumNumbers::SymmetryList const& sl)
{
   using ::CoerceSymmetryList;
   typename ptr_type::lock_type Lock(Data.lock());
   for (iterator I = Lock->begin(); I != Lock->end(); ++I)
   {
      CoerceSymmetryList(I->second, sl);
   }
}

template <typename OperatorT>
OperatorT const&
Block<OperatorT>::operator[](std::string const& s) const
{
   typename DataType::const_iterator I = Data->find(s);
   CHECK(I != Data->end()) << "The block does not contain any operator named " << s;
   return I->second;
}

template <typename OperatorT>
SymmetryList Block<OperatorT>::GetSymmetryList() const
{
   CHECK(!Data->empty());
   // if we're debugging, verify that the symmetry list is the same for all operators
#if !defined(NDEBUG)
   SymmetryList ThisSymmetryList = Data->begin()->second.GetSymmetryList();
   for (const_iterator I = Data->begin(); I != Data->end(); ++I)
   {
      CHECK_EQUAL(I->second.GetSymmetryList(), ThisSymmetryList)(I->first);
   }
#endif
   return Data->begin()->second.GetSymmetryList();
}

template <typename OperatorT>
typename Block<OperatorT>::basis1_type const&
Block<OperatorT>::Basis1() const
{
   CHECK(!Data->empty());
   // if we're debugging, verify that the basis is the same for all operators
#if !defined(NDEBUG)
   basis1_type ThisBasis1 = Data->begin()->second.Basis1();
   for (const_iterator I = Data->begin(); I != Data->end(); ++I)
   {
      CHECK_EQUAL(I->second.Basis1(), ThisBasis1);
   }
#endif
   return Data->begin()->second.Basis1();
}

template <typename OperatorT>
typename Block<OperatorT>::basis2_type const&
Block<OperatorT>::Basis2() const
{
   CHECK(!Data->empty());
   // if we're debugging, verify that the basis is the same for all operators
#if !defined(NDEBUG)
   basis2_type ThisBasis2 = Data->begin()->second.Basis2();
   for (const_iterator I = Data->begin(); I != Data->end(); ++I)
   {
      CHECK_EQUAL(I->second.Basis2(), ThisBasis2);
   }
#endif
   return Data->begin()->second.Basis2();
}

template <typename OperatorT>
PStream::opstream& operator<<(PStream::opstream& out, Block<OperatorT> const& B)
{
   return out << B.Data;
}

template <typename OperatorT>
PStream::ipstream& operator>>(PStream::ipstream& in, Block<OperatorT>& B)
{
   return in >> B.Data;
}
