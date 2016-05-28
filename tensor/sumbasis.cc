// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// tensor/sumbasis.cc
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include <boost/type_traits.hpp>

template <typename BasisT>
inline
int SumBasis<BasisT>::operator()(int BasisNumber, int Subspace) const
{
   DEBUG_RANGE_CHECK_OPEN(BasisNumber, 0, BasisList.size());
   DEBUG_RANGE_CHECK_OPEN(Subspace, 0, BasisList[BasisNumber].size());
   return Mapping[BasisNumber][Subspace];
}

template <>
inline
void SumBasis<SimpleBasis>::AddBasis(SimpleBasis const& B)
{
   PRECONDITION(B.GetSymmetryList() == this->GetSymmetryList());
   BasisList.push_back(B);
   std::vector<int> Map(B.size());
   for (int i = 0; i < B.size(); ++i)
   {
      Map[i] = this->Append(B.qn(i));
   }
   Mapping.push_back(Map);
}

template <>
inline
void SumBasis<VectorBasis>::AddBasis(VectorBasis const& B)
{
   PRECONDITION(B.GetSymmetryList() == this->GetSymmetryList());
   BasisList.push_back(B);
   std::vector<int> Map(B.size());
   for (int i = 0; i < B.size(); ++i)
   {
      Map[i] = this->Append(B.qn(i), B.Dimension(i));
   }
   Mapping.push_back(Map);
}


template <typename BasisT>
SumBasis<BasisT>::SumBasis(SymmetryList const& SList)
   : BasisType(SList)
{
}

template <typename BasisT>
SumBasis<BasisT>::SumBasis(BasisType const& B)
   : BasisType(B.GetSymmetryList())
{
   this->AddBasis(B);
}
   
template <typename BasisT>
SumBasis<BasisT>::SumBasis(BasisType const& B1, BasisType const& B2)
   : BasisType(B1.GetSymmetryList())
{
   this->AddBasis(B1);
   this->AddBasis(B2);
}


template <typename C, typename Basis1T, typename Basis2T>
IrredOperator<C, Basis1T, Basis2T>
tensor_sum(IrredOperator<C, Basis1T, Basis2T> const& x,
	   IrredOperator<C, Basis1T, Basis2T> const& y, 
           SumBasis<Basis1T> const& B1, 
           SumBasis<Basis2T> const& B2)
{
   DEBUG_PRECONDITION(x.TransformsAs() == y.TransformsAs());
   DEBUG_PRECONDITION(B1.NumBasis() == 2);
   DEBUG_PRECONDITION(B2.NumBasis() == 2);
   DEBUG_PRECONDITION(B1.Basis(0) == x.Basis1());
   DEBUG_PRECONDITION(B1.Basis(1) == y.Basis1());
   DEBUG_PRECONDITION(B2.Basis(0) == x.Basis2());
   DEBUG_PRECONDITION(B2.Basis(1) == y.Basis2());

   typedef typename IrredOperator<C, Basis1T, Basis2T>::const_iterator1 const_iterator1;
   IrredOperator<C, Basis1T, Basis2T> Result(B1, B2, x.TransformsAs());
   int i = 0;
   for (const_iterator1 I = x.begin1(); I != x.end1(); ++I, ++i)
   {
      for (typename const_iterator1::const_iterator J = I.begin(); J != I.end(); ++J)
      {
	 Result(B1(0,i), B2(0,J.index())) = J.data();
      }
   }
   i = 0;
   for (const_iterator1 I = y.begin1(); I != y.end1(); ++I, ++i)
   {
      for (typename const_iterator1::const_iterator J = I.begin(); J != I.end(); ++J)
      {
	 Result(B1(1,i), B2(1,J.index())) = J.data();
      }
   }
   return Result;
}

template <typename C, typename Basis1T, typename Basis2T>
inline
IrredOperator<C, Basis1T, Basis2T>
tensor_sum(IrredOperator<C, Basis1T, Basis2T> const& x,
	   IrredOperator<C, Basis1T, Basis2T> const& y)
{
   return tensor_sum(x, y, 
                     SumBasis<Basis1T>(x.Basis1(), y.Basis1()), 
                     SumBasis<Basis2T>(x.Basis2(), y.Basis2()));
}

template <typename C, typename Basis1T, typename Basis2T>
IrredOperator<C, Basis1T, Basis2T>
tensor_row_sum(IrredOperator<C, Basis1T, Basis2T> const& x,
	       IrredOperator<C, Basis1T, Basis2T> const& y, SumBasis<Basis2T> const& B2)
{
   DEBUG_PRECONDITION(x.TransformsAs() == y.TransformsAs());
   DEBUG_PRECONDITION(B2.NumBasis() == 2);
   DEBUG_PRECONDITION(x.Basis1() == y.Basis1());
   DEBUG_PRECONDITION(x.Basis2() == B2.Basis(0));
   DEBUG_PRECONDITION(y.Basis2() == B2.Basis(1));

   typedef typename IrredOperator<C, Basis1T, Basis2T>::const_iterator1 const_iterator1;
   IrredOperator<C, Basis1T, Basis2T> Result(x.Basis1(), B2, x.TransformsAs());
   int i = 0;
   for (const_iterator1 I = x.begin1(); I != x.end1(); ++I, ++i)
   {
      for (typename const_iterator1::const_iterator J = I.begin(); J != I.end(); ++J)
      {
	 Result(i, B2(0,J.index())) = J.data();
      }
   }
   i = 0;
   for (const_iterator1 I = y.begin1(); I != y.end1(); ++I, ++i)
   {
      for (typename const_iterator1::const_iterator J = I.begin(); J != I.end(); ++J)
      {
	 Result(i, B2(1,J.index())) = J.data();
      }
   }
   return Result;
}

template <typename C, typename Basis1T, typename Basis2T>
inline
IrredOperator<C, Basis1T, Basis2T>
tensor_row_sum(IrredOperator<C, Basis1T, Basis2T> const& x,
	       IrredOperator<C, Basis1T, Basis2T> const& y)
{
   return tensor_row_sum(x, y, SumBasis<Basis2T>(x.Basis2(), y.Basis2()));
}

template <typename C, typename Basis1T, typename Basis2T>
IrredOperator<C, Basis1T, Basis2T>
tensor_col_sum(IrredOperator<C, Basis1T, Basis2T> const& x,
	       IrredOperator<C, Basis1T, Basis2T> const& y, SumBasis<Basis1T> const& B1)
{
   DEBUG_PRECONDITION(x.TransformsAs() == y.TransformsAs());
   DEBUG_PRECONDITION(B1.NumBasis() == 2);
   DEBUG_PRECONDITION(x.Basis2() == y.Basis2());
   DEBUG_PRECONDITION(x.Basis1() == B1.Basis(0));
   DEBUG_PRECONDITION(y.Basis1() == B1.Basis(1));

   typedef typename IrredOperator<C, Basis1T, Basis2T>::const_iterator1 const_iterator1;
   IrredOperator<C, Basis1T, Basis2T> Result(B1, x.Basis2(), x.TransformsAs());
   int i = 0;
   for (const_iterator1 I = x.begin1(); I != x.end1(); ++I, ++i)
   {
      for (typename const_iterator1::const_iterator J = I.begin(); J != I.end(); ++J)
      {
	 Result(B1(0,i), J.index()) = J.data();
      }
   }
   i = 0;
   for (const_iterator1 I = y.begin1(); I != y.end1(); ++I, ++i)
   {
      for (typename const_iterator1::const_iterator J = I.begin(); J != I.end(); ++J)
      {
	 Result(B1(1,i), J.index()) = J.data();
      }
   }
   return Result;
}

template <typename C, typename Basis1T, typename Basis2T>
inline
IrredOperator<C, Basis1T, Basis2T>
tensor_col_sum(IrredOperator<C, Basis1T, Basis2T> const& x,
	       IrredOperator<C, Basis1T, Basis2T> const& y)
{
   return tensor_col_sum(x, y, SumBasis<Basis1T>(x.Basis1(), y.Basis1()));
}




template <typename FwdIter, typename Basis1T, typename Basis2T>
typename std::iterator_traits<FwdIter>::value_type
tensor_accumulate(FwdIter first, FwdIter last,
           SumBasis<Basis1T> const& B1, 
           SumBasis<Basis2T> const& B2)
{
   typedef typename std::iterator_traits<FwdIter>::value_type ResultType;
   BOOST_STATIC_ASSERT((boost::is_same<typename ResultType::Basis1Type, Basis1T>::value));
   BOOST_STATIC_ASSERT((boost::is_same<typename ResultType::Basis2Type, Basis2T>::value));

   FwdIter F = first;
   int b = 0;
   while (F != last)
   {
      CHECK_EQUAL(F->TransformsAs(), first->TransformsAs());
      CHECK_EQUAL(F->Basis2(), B2.Basis(b));
      CHECK_EQUAL(F->Basis1(), B1.Basis(b));
      ++b;
      ++F;
   }
   CHECK_EQUAL(B1.NumBasis(), b);
   CHECK_EQUAL(B2.NumBasis(), b);

   typedef typename ResultType::const_iterator1 const_iterator1;
   ResultType Result(B1, B2, first->TransformsAs());
   int f = 0;
   for (F = first; F != last; ++F, ++f)
   {
      int i = 0;
      ResultType const Element = *F;
      for (const_iterator1 I = Element.begin1(); I != Element.end1(); ++I, ++i)
      {
	 for (typename const_iterator1::const_iterator J = I.begin(); J != I.end(); ++J)
	 {
	    Result(B1(f,i), B2(f,J.index())) = J.data();
	 }
      }
   }
   return Result;
}

template <typename FwdIter, typename Basis2T>
typename std::iterator_traits<FwdIter>::value_type
tensor_row_accumulate(FwdIter first, FwdIter last, SumBasis<Basis2T> const& B2)
{
   typedef typename std::iterator_traits<FwdIter>::value_type ResultType;
   BOOST_STATIC_ASSERT((boost::is_same<typename ResultType::Basis2Type, Basis2T>::value));

   CHECK_EQUAL(first->Basis2(), B2.Basis(0));
   FwdIter F = first;
   ++F;
   int b = 1;
   while (F != last)
   {
      CHECK_EQUAL(F->TransformsAs(), first->TransformsAs());
      CHECK_EQUAL(F->Basis1(), first->Basis1());
      CHECK_EQUAL(F->Basis2(), B2.Basis(b));
      ++b;
      ++F;
   }
   CHECK_EQUAL(B2.NumBasis(), b);

   typedef typename ResultType::const_iterator1 const_iterator1;
   ResultType Result(first->Basis1(), B2, first->TransformsAs());
   int f = 0;
   for (F = first; F != last; ++F, ++f)
   {
      int i = 0;
      ResultType const Element = *F;
      for (const_iterator1 I = Element.begin1(); I != Element.end1(); ++I, ++i)
      {
	 for (typename const_iterator1::const_iterator J = I.begin(); J != I.end(); ++J)
	 {
	    Result(i, B2(f,J.index())) = J.data();
	 }
      }
   }
   return Result;
}

template <typename FwdIter, typename Basis1T>
typename std::iterator_traits<FwdIter>::value_type
tensor_col_accumulate(FwdIter first, FwdIter last, SumBasis<Basis1T> const& B1)
{
   typedef typename std::iterator_traits<FwdIter>::value_type ResultType;
   BOOST_STATIC_ASSERT((boost::is_same<typename ResultType::Basis1Type, Basis1T>::value));

   if (first == last) return ResultType();

   CHECK_EQUAL(first->Basis1(), B1.Basis(0));
   FwdIter F = first;
   ++F;
   int b = 1;
   while (F != last)
   {
      CHECK_EQUAL(F->TransformsAs(), first->TransformsAs());
      CHECK_EQUAL(F->Basis2(), first->Basis2());
      CHECK_EQUAL(F->Basis1(), B1.Basis(b));
      ++b;
      ++F;
   }
   CHECK_EQUAL(B1.NumBasis(), b);

   typedef typename ResultType::const_iterator1 const_iterator1;
   ResultType Result(B1, first->Basis2(), first->TransformsAs());
   int f = 0;
   for (F = first; F != last; ++F, ++f)
   {
      int i = 0;
      ResultType const Element = *F;
      for (const_iterator1 I = Element.begin1(); I != Element.end1(); ++I, ++i)
      {
	 for (typename const_iterator1::const_iterator J = I.begin(); J != I.end(); ++J)
	 {
	    Result(B1(f,i), J.index()) = J.data();
	 }
      }
   }
   return Result;
}
