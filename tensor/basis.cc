// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// tensor/basis.cc
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

#include <iomanip>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <boost/iterator/transform_iterator.hpp>

namespace Tensor
{

//
// BasisList
//

inline
PStream::opstream& operator<<(PStream::opstream& out, BasisList const& B)
{
   return out << B.S_ << B.Q_;
}

inline
PStream::ipstream& operator>>(PStream::ipstream& in, BasisList& B)
{
   return in >> B.S_ >> B.Q_;
}

template <typename FwdIter>
BasisList::BasisList(FwdIter first, FwdIter last)
   : S_(first == last ? QuantumNumbers::SymmetryList() : first->GetSymmetryList()), Q_(first, last)
{
   DEBUG_CHECK(first != last)("The list must be non-empty");
}

template <typename FwdIter>
BasisList::BasisList(QuantumNumbers::SymmetryList const& S, FwdIter first, FwdIter last)
   : S_(S), Q_(first, last)
{
}
inline
BasisList adjoint(BasisList const& b)
{
   typedef Adjoint<QuantumNumber> Adj;
   typedef boost::transform_iterator<Adj, BasisList::const_iterator> AdjointIter;
   return BasisList(AdjointIter(b.begin(), Adj()), AdjointIter(b.end(), Adj()));
}

inline
BasisList
make_vacuum_basis(SymmetryList const& S)
{
   BasisList Result(S);
   Result.push_back(QuantumNumber(S));
   return Result;
}

inline
BasisList
make_single_basis(QuantumNumbers::QuantumNumber const& q)
{
   BasisList Result(q.GetSymmetryList());
   Result.push_back(q);
   return Result;
}

//
// VectorBasis
//

inline
PStream::opstream& operator<<(PStream::opstream& out, VectorBasis const& B)
{
   return out << B.Basis_ << B.Dimension_;
}

inline
PStream::ipstream& operator>>(PStream::ipstream& in, VectorBasis& B)
{
   return in >> B.Basis_ >> B.Dimension_;
}

inline
VectorBasis::VectorBasis(SymmetryList const& sl)
   : Basis_(sl)
{
}

inline
VectorBasis::VectorBasis(BasisList const& Basis)
   : Basis_(Basis), Dimension_(Basis.size(), 1)
{
}

template <typename FwdIter>
VectorBasis::VectorBasis(BasisList const& Basis, FwdIter first, FwdIter last)
   : Basis_(Basis), Dimension_(first, last)
{
}


template <typename U>
struct select1st_from_iter
{
   typedef typename std::iterator_traits<U>::value_type T;
   typedef T const& argument_type;
   typedef typename T::first_type result_type;

   result_type operator()(argument_type x) const
   {
      return x.first;
   }
};

template <typename U>
struct select2nd_from_iter
{
   typedef typename std::iterator_traits<U>::value_type T;
   typedef T const& argument_type;
   typedef typename T::second_type result_type;

   result_type operator()(argument_type x) const
   {
      return x.second;
   }
};

template <typename FwdIter>
VectorBasis::VectorBasis(FwdIter first, FwdIter last)
   : Basis_(boost::make_transform_iterator(first, select1st_from_iter<FwdIter>()),
            boost::make_transform_iterator(last, select1st_from_iter<FwdIter>())),
     Dimension_(boost::make_transform_iterator(first, select2nd_from_iter<FwdIter>()),
                boost::make_transform_iterator(last, select2nd_from_iter<FwdIter>()))
{
}

inline
void VectorBasis::push_back(QuantumNumber const& q, int Dimension)
{
   Basis_.push_back(q);
   Dimension_.push_back(Dimension);
}

inline
int VectorBasis::total_dimension() const
{
   return std::accumulate(Dimension_.begin(), Dimension_.end(), 0);
}

inline
VectorBasis adjoint(VectorBasis const& b)
{
   return VectorBasis(adjoint(b.Basis_), b.Dimension_);
}

} // namespace Tensor
