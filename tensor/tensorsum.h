// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// tensor/tensorsum.h
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

#if !defined(TENSORSUM_H_KJSAHEWHOIUHIUFHAWELO9H)
#define TENSORSUM_H_KJSAHEWHOIUHIUFHAWELO9H

#include "tensor.h"

namespace Tensor
{

template <typename T>
class SumBasis;

template <>
class SumBasis<BasisList> : public BasisList
{
   public:
      typedef int target_type;
      // in the opposite direction, we map a subspace into a pair of original subspaces
      typedef std::pair<int, int> source_type;

   private:
      typedef std::list<target_type> TargetListType;
      typedef LinearAlgebra::Vector<source_type>    SourceListType;

   public:
      typedef TargetListType::const_iterator const_iterator;

      SumBasis() {}

      explicit SumBasis(SymmetryList const& SList);
      explicit SumBasis(BasisList const& B);
      SumBasis(BasisList const& Basis1_, BasisList const& Basis2_);

      template <typename FwdIter>
      SumBasis(FwdIter first, FwdIter last);

      int NumBasis() const { return BList_.size(); }

      void AddBasis(BasisList const& B);

      int operator()(int BasisNumber, int Subspace) const;

      BasisList const& Basis(int BasisNumber) const { return BList_[BasisNumber]; }

      BasisList const& Basis() const { return *this; }

   private:
      BasisList& base() { return *this; }

      std::vector<BasisList> BList_;
      std::vector<std::vector<int> > Mapping_;

      // we probably don't need the reverse mapping
      // std::vector<source_type> ReverseMapping;

   friend PStream::opstream& operator<<(PStream::opstream& out, SumBasis<BasisList> const& B);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, SumBasis<BasisList>& B);
};

template <>
class SumBasis<VectorBasis>
{
   public:
      typedef VectorBasis::value_type value_type;

      SumBasis() {}

      explicit SumBasis(SymmetryList const& SList);
      explicit SumBasis(VectorBasis const& B);

      SumBasis(VectorBasis const& Basis1, VectorBasis const& Basis2);

      template <typename FwdIter>
      SumBasis(FwdIter first, FwdIter last);

      int NumBasis() const { return BList_.size(); }

      void AddBasis(VectorBasis const& B);

      int operator()(int BasisNumber, int Subspace) const;

      operator VectorBasis const&() const { return This_; }

      VectorBasis const& Basis() const { return This_; }

      VectorBasis const& Basis(int bn) const { return BList_[bn]; }

      SymmetryList const& GetSymmetryList() const { return This_.GetSymmetryList(); }

      value_type const& operator[](int s) const { return This_[s]; }
      int dim(int s) const { return This_.dim(s); }

      int total_dimension() const { return This_.total_dimension(); }

      int total_degree() const { return This_.total_degree(); }

   private:
      VectorBasis This_;

      std::vector<VectorBasis> BList_;
      std::vector<std::vector<int> > Mapping_;

   friend PStream::opstream& operator<<(PStream::opstream& out, SumBasis<VectorBasis> const& B);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, SumBasis<VectorBasis>& B);
};

//
// tensor_sum
//
// returns the direct sum of irreducible tensor operators x and y.
// Precondition: x.TransformsAs() == y.TransformsAs()
//            && x.Basis1() == B1.Basis(0)
//            && y.Basis1() == B1.Basis(1)
//            && x.Basis2() == B2.Basis(0)
//            && y.Basis2() == B2.Basis(1)

template <typename T, typename B1, typename B2, typename S>
IrredTensor<T, B1, B2, S>
tensor_sum(IrredTensor<T, B1, B2, S> const& x, IrredTensor<T, B1, B2, S> const& y,
           SumBasis<B1> const& b1, SumBasis<B2> const& b2);

template <typename T, typename B, typename S>
inline
IrredTensor<T, B, B, S>
tensor_sum(IrredTensor<T, B, B, S> const& x, IrredTensor<T, B, B, S> const& y,
           SumBasis<B> const& b)
{
   return tensor_sum(x,y,b,b);
}

template <typename T, typename B1, typename B2, typename S>
IrredTensor<T, B1, B2, S>
tensor_sum(IrredTensor<T, B1, B2, S> const& x, IrredTensor<T, B1, B2, S> const& y);

// constructs the 'horizontal concatenation' of x and y,
// ie. a row-vector (x,y)
template <typename T, typename B1, typename B2, typename S>
IrredTensor<T, B1, B2, S>
tensor_row_sum(IrredTensor<T, B1, B2, S> const& x,
               IrredTensor<T, B1, B2, S> const& y,
               SumBasis<B2> const& b2);

template <typename T, typename B1, typename B2, typename S>
inline
IrredTensor<T, B1, B2, S>
tensor_row_sum(IrredTensor<T, B1, B2, S> const& x, IrredTensor<T, B1, B2, S> const& y);

// constructs the 'vertical concatenation' of x and y,
// ie. a column-vector (x)
//                     (y)
template <typename T, typename B1, typename B2, typename S>
IrredTensor<T, B1, B2, S>
tensor_col_sum(IrredTensor<T, B1, B2, S> const& x, IrredTensor<T, B1, B2, S> const& y,
               SumBasis<B1> const& b1);

template <typename T, typename B1, typename B2, typename S>
IrredTensor<T, B1, B2, S>
tensor_col_sum(IrredTensor<T, B1, B2, S> const& x, IrredTensor<T, B1, B2, S> const& y);

template <typename FwdIter, typename B1, typename B2>
typename std::iterator_traits<FwdIter>::value_type
tensor_accumulate(FwdIter first, FwdIter last,
           SumBasis<B1> const& b1,
           SumBasis<B2> const& b2);

template <typename FwdIter, typename B2>
typename std::iterator_traits<FwdIter>::value_type
tensor_row_accumulate(FwdIter first, FwdIter last, SumBasis<B2> const& b2);

template <typename FwdIter, typename B1>
typename std::iterator_traits<FwdIter>::value_type
tensor_col_accumulate(FwdIter first, FwdIter last, SumBasis<B1> const& b1);

} // namespace Tensor

#include "tensorsum.cc"

#endif
