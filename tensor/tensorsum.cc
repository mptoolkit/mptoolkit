// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// tensor/tensorsum.cc
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

// SumBasis<BasisList>

inline
SumBasis<BasisList>::SumBasis(SymmetryList const& SList)
   : BasisList(SList)
{
}

inline
SumBasis<BasisList>::SumBasis(BasisList const& B)
   : BasisList(B.GetSymmetryList())
{
   this->AddBasis(B);
}

inline
SumBasis<BasisList>::SumBasis(BasisList const& B1, BasisList const& B2)
   : BasisList(B1.GetSymmetryList())
{
   this->AddBasis(B1);
   this->AddBasis(B2);
}

template <typename FwdIter>
SumBasis<BasisList>::SumBasis(FwdIter first, FwdIter last)
   : BasisList(first->GetSymmetryList())
{
   while (first != last)
   {
      this->AddBasis(*first);
      ++first;
   }
}

inline
int SumBasis<BasisList>::operator()(int BasisNumber, int Subspace) const
{
   DEBUG_RANGE_CHECK_OPEN(BasisNumber, 0, int(BList_.size()));
   DEBUG_RANGE_CHECK_OPEN(Subspace, 0, int(BList_[BasisNumber].size()));
   return Mapping_[BasisNumber][Subspace];
}

// SumBasis<VectorBasis>

inline
SumBasis<VectorBasis>::SumBasis(SymmetryList const& SList)
   : This_(SList)
{
}

inline
SumBasis<VectorBasis>::SumBasis(VectorBasis const& B)
   : This_(B.GetSymmetryList())
{
   this->AddBasis(B);
}

inline
SumBasis<VectorBasis>::SumBasis(VectorBasis const& B1, VectorBasis const& B2)
   : This_(B1.GetSymmetryList())
{
   this->AddBasis(B1);
   this->AddBasis(B2);
}

template <typename FwdIter>
SumBasis<VectorBasis>::SumBasis(FwdIter first, FwdIter last)
   : This_(first->GetSymmetryList())
{
   while (first != last)
   {
      this->AddBasis(*first);
      ++first;
   }
}

inline
int SumBasis<VectorBasis>::operator()(int BasisNumber, int Subspace) const
{
   DEBUG_RANGE_CHECK_OPEN(BasisNumber, 0, int(BList_.size()));
   DEBUG_RANGE_CHECK_OPEN(Subspace, 0, int(BList_[BasisNumber].size()));
   return Mapping_[BasisNumber][Subspace];
}

// tensor_sum

template <typename T, typename B1, typename B2, typename S>
inline
IrredTensor<T, B1, B2, S>
tensor_sum(IrredTensor<T, B1, B2, S> const& x, IrredTensor<T, B1, B2, S> const& y)
{
   return tensor_sum(x, y,
                     SumBasis<B1>(x.Basis1(), y.Basis1()),
                     SumBasis<B2>(x.Basis2(), y.Basis2()));
}


template <typename T, typename B1, typename B2, typename S>
inline
IrredTensor<T, B1, B2, S>
tensor_row_sum(IrredTensor<T, B1, B2, S> const& x, IrredTensor<T, B1, B2, S> const& y)
{
   return tensor_row_sum(x, y, SumBasis<B2>(x.Basis2(), y.Basis2()));
}

template <typename T, typename B1, typename B2, typename S>
inline
IrredTensor<T, B1, B2, S>
tensor_col_sum(IrredTensor<T, B1, B2, S> const& x, IrredTensor<T, B1, B2, S> const& y)
{
   return tensor_col_sum(x, y, SumBasis<B1>(x.Basis1(), y.Basis1()));
}

template <typename T, typename B1, typename B2, typename S>
IrredTensor<T, B1, B2, S>
tensor_sum(IrredTensor<T, B1, B2, S> const& x, IrredTensor<T, B1, B2, S> const& y,
           SumBasis<B1> const& b1, SumBasis<B2> const& b2)
{
   DEBUG_PRECONDITION(x.TransformsAs() == y.TransformsAs());
   DEBUG_PRECONDITION(b1.NumBasis() == 2);
   DEBUG_PRECONDITION(b2.NumBasis() == 2);
   DEBUG_PRECONDITION(b1.Basis(0) == x.Basis1());
   DEBUG_PRECONDITION(b1.Basis(1) == y.Basis1());
   DEBUG_PRECONDITION(b2.Basis(0) == x.Basis2());
   DEBUG_PRECONDITION(b2.Basis(1) == y.Basis2());

   typedef typename const_iterator<IrredTensor<T, S> >::type iterator;
   typedef typename const_inner_iterator<IrredTensor<T, S> >::type inner_iterator;
   IrredTensor<T, B1, B2, S> Result(b1, b2, x.TransformsAs());
   for (iterator I = iterate(x); I; ++I)
   {
      for (inner_iterator J = iterate(I); J; ++J)
      {
         Result(b1(0,J.index1()), b2(0,J.index2())) = *J;
      }
   }
   for (iterator I = iterate(y); I; ++I)
   {
      for (inner_iterator J = iterate(I); J; ++J)
      {
         Result(b1(1,J.index1()), b2(1,J.index2())) = *J;
      }
   }
   return Result;
}

template <typename T, typename B1, typename B2, typename S>
IrredTensor<T, B1, B2, S>
tensor_row_sum(IrredTensor<T, B1, B2, S> const& x,
               IrredTensor<T, B1, B2, S> const& y, SumBasis<B2> const& b2)
{
   DEBUG_PRECONDITION_EQUAL(x.TransformsAs(), y.TransformsAs());
   DEBUG_PRECONDITION_EQUAL(b2.NumBasis(), 2);
   DEBUG_PRECONDITION_EQUAL(x.Basis1(), y.Basis1());
   DEBUG_PRECONDITION_EQUAL(x.Basis2(), b2.Basis(0));
   DEBUG_PRECONDITION_EQUAL(y.Basis2(), b2.Basis(1));

   typedef typename const_iterator<IrredTensor<T, B1, B2, S> >::type iterator;
   typedef typename const_inner_iterator<IrredTensor<T, B1, B2, S> >::type inner_iterator;
   IrredTensor<T, B1, B2, S> Result(x.Basis1(), b2, x.TransformsAs());
   for (iterator I = iterate(x); I; ++I)
   {
      for (inner_iterator J = iterate(I); J; ++J)
      {
         Result(J.index1(), b2(0,J.index2())) = *J;
      }
   }
   for (iterator I = iterate(y); I; ++I)
   {
      for (inner_iterator J = iterate(I); J; ++J)
      {
         Result(J.index1(), b2(1,J.index2())) = *J;
      }
   }
   return Result;
}

template <typename T, typename B1, typename B2, typename S>
IrredTensor<T, B1, B2, S>
tensor_col_sum(IrredTensor<T, B1, B2, S> const& x,
               IrredTensor<T, B1, B2, S> const& y,
               SumBasis<B1> const& b1)
{
   DEBUG_PRECONDITION_EQUAL(x.TransformsAs(), y.TransformsAs());
   DEBUG_PRECONDITION_EQUAL(b1.NumBasis(), 2);
   DEBUG_PRECONDITION_EQUAL(x.Basis2(), y.Basis2());
   DEBUG_PRECONDITION_EQUAL(x.Basis1(), b1.Basis(0));
   DEBUG_PRECONDITION_EQUAL(y.Basis1(), b1.Basis(1));

   typedef typename const_iterator<IrredTensor<T, B1, B2, S> >::type iterator;
   typedef typename const_inner_iterator<IrredTensor<T, B1, B2, S> >::type inner_iterator;
   IrredTensor<T, B1, B2, S> Result(b1, x.Basis2(), x.TransformsAs());
   for (iterator I = iterate(x); I; ++I)
   {
      for (inner_iterator J = iterate(I); J; ++J)
      {
         Result(b1(0,J.index1()), J.index2()) = *J;
      }
   }
   for (iterator I = iterate(y); I; ++I)
   {
      for (inner_iterator J = iterate(I); J; ++J)
      {
         Result(b1(1,J.index1()), J.index2()) = *J;
      }
   }
   return Result;
}

template <typename FwdIter, typename B1, typename B2>
typename std::iterator_traits<FwdIter>::value_type
tensor_accumulate(FwdIter first, FwdIter last,
                  SumBasis<B1> const& b1,
                  SumBasis<B2> const& b2)
{
   typedef typename std::iterator_traits<FwdIter>::value_type ResultType;

   FwdIter F = first;
   int b = 0;
   while (F != last)
   {
      CHECK_EQUAL(F->TransformsAs(), first->TransformsAs());
      CHECK_EQUAL(F->Basis2(), b2.Basis(b));
      CHECK_EQUAL(F->Basis1(), b1.Basis(b));
      ++b;
      ++F;
   }
   CHECK_EQUAL(b1.NumBasis(), b);
   CHECK_EQUAL(b2.NumBasis(), b);

   typedef typename const_iterator<ResultType>::type iterator;
   typedef typename const_inner_iterator<ResultType>::type inner_iterator;
   ResultType Result(b1, b2, first->TransformsAs());
   int f = 0;
   for (F = first; F != last; ++F, ++f)
   {
      for (iterator I = iterate(*F); I; ++I)
      {
         for (inner_iterator J = iterate(I); J; ++J)
         {
            Result(b1(f,J.index1()), b2(f,J.index2())) = *J;
         }
      }
   }
   return Result;
}

template <typename FwdIter, typename B2>
typename std::iterator_traits<FwdIter>::value_type
tensor_row_accumulate(FwdIter first, FwdIter last, SumBasis<B2> const& b2)
{
   typedef typename std::iterator_traits<FwdIter>::value_type ResultType;

   CHECK_EQUAL(first->Basis2(), b2.Basis(0));
   FwdIter F = first;
   ++F;
   int b = 1;
   while (F != last)
   {
      CHECK_EQUAL(F->TransformsAs(), first->TransformsAs());
      CHECK_EQUAL(F->Basis1(), first->Basis1());
      CHECK_EQUAL(F->Basis2(), b2.Basis(b));
      ++b;
      ++F;
   }
   CHECK_EQUAL(b2.NumBasis(), b);

   typedef typename const_iterator<ResultType>::type iterator;
   typedef typename const_inner_iterator<ResultType>::type inner_iterator;
   ResultType Result(first->Basis1(), b2, first->TransformsAs());
   int f = 0;
   for (F = first; F != last; ++F, ++f)
   {
      for (iterator I = iterate(*F); I; ++I)
      {
         for (inner_iterator J = iterate(I); J; ++J)
         {
            Result(J.index1(), b2(f,J.index2())) = *J;
         }
      }
   }
   return Result;
}

template <typename FwdIter, typename B1>
typename std::iterator_traits<FwdIter>::value_type
tensor_col_accumulate(FwdIter first, FwdIter last, SumBasis<B1> const& b1)
{
   typedef typename std::iterator_traits<FwdIter>::value_type ResultType;

   if (first == last) return ResultType();

   CHECK_EQUAL(first->Basis1(), b1.Basis(0));
   FwdIter F = first;
   ++F;
   int b = 1;
   while (F != last)
   {
      CHECK_EQUAL(F->TransformsAs(), first->TransformsAs());
      CHECK_EQUAL(F->Basis2(), first->Basis2());
      CHECK_EQUAL(F->Basis1(), b1.Basis(b));
      ++b;
      ++F;
   }
   CHECK_EQUAL(b1.NumBasis(), b);

   typedef typename const_iterator<ResultType>::type iterator;
   typedef typename const_inner_iterator<ResultType>::type inner_iterator;
   ResultType Result(b1, first->Basis2(), first->TransformsAs());
   int f = 0;
   for (F = first; F != last; ++F, ++f)
   {
      for (iterator I = iterate(*F); I; ++I)
      {
         for (inner_iterator J = iterate(I); J; ++J)
         {
            Result(b1(f,J.index1()), J.index2()) = *J;
         }
      }
   }
   return Result;
}

} // namespace Tensor
