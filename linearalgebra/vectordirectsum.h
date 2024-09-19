// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/vectordirectsum.h
//
// Copyright (C) 2005-2016 Ian McCulloch <ian@qusim.net>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

/*
  vectordirectsum.h

  implementation of VectorDirectSum.

  Created 2005-03-16 Ian McCulloch
*/

#if !defined(VECTORDIRECTSUM_H_DSRIGHTUY389Y89Y89PYP)
#define VECTORDIRECTSUM_H_DSRIGHTUY389Y89Y89PYP

#include "vectoroperationsbase.h"

namespace LinearAlgebra
{

// TODO: this is a dumb implementation, should return a proxy.

// TODO: these should use enable_if on is_defined<Addition<Sv, Tv> >

template <typename S, typename T, typename Sv, typename Si,
          typename Tv, typename Ti>
struct VectorDirectSum<S, T, LOCAL_VECTOR(Sv, Si),  LOCAL_VECTOR(Tv, Ti)>
{
   typedef typename vector_abstract_or<typename abstract_interface<S>::type,
                                       typename abstract_interface<T>::type>::type ai;

   typedef typename make_value<typename Addition<Sv, Tv>::result_type>::type value_type;
   typedef typename make_vector_from_abstract<value_type, ai>::type result_type;
   typedef S const& first_argument_type;
   typedef T const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      result_type R(size(x) + size(y));
      zero_all(R);
      typename const_iterator<S>::type I = iterate(x);
      while (I)
      {
         set_element(R, I.index(), *I);
         ++I;
      }
      size_type xSize = size(x);
      typename const_iterator<T>::type J = iterate(y);
      while (J)
      {
         set_element(R, xSize + J.index(), *J);
         ++J;
      }
      return R;
   }
};

template <typename S, typename T, typename Sv, typename Si,
          typename Tv, typename Ti>
struct VectorDirectSum<S, T, DENSE_VECTOR(Sv, Si),  DENSE_VECTOR(Tv, Ti)>
{
   typedef typename vector_abstract_or<typename abstract_interface<S>::type,
                                       typename abstract_interface<T>::type>::type ai;

   typedef typename make_value<typename Addition<Sv, Tv>::result_type>::type value_type;
   typedef typename make_vector_from_abstract<value_type, ai>::type result_type;
   typedef S const& first_argument_type;
   typedef T const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      result_type R(size(x) + size(y));
      typename iterator<result_type>::type K = iterate(R);
      typename const_iterator<S>::type I = iterate(x);
      while (I)
      {
         *K = *I;
         ++I;
         ++K;
      }
      typename const_iterator<T>::type J = iterate(y);
      while (J)
      {
         *K = *J;
         ++J;
         ++K;
      }
      return R;
   }
};

// specialization doesn't require Addition to be defined
template <typename S, typename T, typename Val, typename Si,
          typename Ti>
struct VectorDirectSum<S, T, LOCAL_VECTOR(Val, Si),  LOCAL_VECTOR(Val, Ti)>
{
   typedef typename vector_abstract_or<typename abstract_interface<S>::type,
                                       typename abstract_interface<T>::type>::type ai;
   typedef typename make_vector_from_abstract<Val, ai>::type result_type;
   typedef S const& first_argument_type;
   typedef T const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      result_type R(size(x) + size(y));
      zero_all(R);
      typename const_iterator<S>::type I = iterate(x);
      while (I)
      {
         set_element(R, I.index(), *I);
         ++I;
      }
      size_type xSize = size(x);
      typename const_iterator<T>::type J = iterate(y);
      while (J)
      {
         set_element(R, xSize + J.index(), *J);
         ++J;
      }
      return R;
   }
};

template <typename S, typename T, typename Val, typename Si,
          typename Ti>
struct VectorDirectSum<S, T, DENSE_VECTOR(Val, Si),  DENSE_VECTOR(Val, Ti)>
{
   typedef typename vector_abstract_or<typename abstract_interface<S>::type,
                                       typename abstract_interface<T>::type>::type ai;
   typedef typename make_vector_from_abstract<Val, ai>::type result_type;
   typedef S const& first_argument_type;
   typedef T const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      result_type R(size(x) + size(y));
      typename iterator<result_type>::type K = iterate(R);
      typename const_iterator<S>::type I = iterate(x);
      while (I)
      {
         *K = *I;
         ++I;
         ++K;
      }
      typename const_iterator<T>::type J = iterate(y);
      while (J)
      {
         *K = *J;
         ++J;
         ++K;
      }
      return R;
   }
};

// TODO: this should use enable_if on is_defined<Addition<Sv, Tv> >
template <typename S, typename T, typename Sv, typename Si, typename Tv>
struct DirectSumInterface<S, T, LOCAL_VECTOR(Sv, Si),  AnyScalar<Tv> >
{
   typedef typename make_value<typename Addition<Sv, Tv>::result_type>::type value_type;

   typedef typename make_vector_from_abstract<value_type,
      typename abstract_interface<S>::type>::type result_type;

   typedef S const& first_argument_type;
   typedef T const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      result_type R(size(x) + 1);
      zero_all(R);
      typename const_iterator<S>::type I = iterate(x);
      while (I)
      {
         set_element(R, I.index(), *I);
         ++I;
      }
      set_element(R, size(x), y);
      return R;
   }
};

// specialization doesn't require Addition to be defined
template <typename S, typename T, typename Val, typename Si>
struct DirectSumInterface<S, T, LOCAL_VECTOR(Val, Si),  AnyScalar<Val> >
{
   typedef typename make_value<S>::type result_type;
   typedef S const& first_argument_type;
   typedef T const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      result_type R(size(x) + 1);
      zero_all(R);
      typename const_iterator<S>::type I = iterate(x);
      while (I)
      {
         set_element(R, I.index(), *I);
         ++I;
      }
      set_element(R, size(x), y);
      return R;
   }
};

template <typename S, typename T, typename Val, typename Si>
struct DirectSumInterface<S, T, DENSE_VECTOR(Val, Si),  AnyScalar<Val> >
{
   typedef typename make_value<S>::type result_type;
   typedef S const& first_argument_type;
   typedef T const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      result_type R(size(x) + 1);
      typename const_iterator<S>::type I = iterate(x);
      while (I)
      {
         set_element(R, I.index(), *I);
         ++I;
      }
      set_element(R, size(x), y);
      return R;
   }
};

template <typename S, typename T, typename Sv, typename Tv, typename Ti>
struct DirectSumInterface<S, T, AnyScalar<Sv>, LOCAL_VECTOR(Tv, Ti)>
{
   typedef typename make_value<typename Addition<Sv, Tv>::result_type>::type value_type;

   typedef typename make_vector_from_abstract<value_type,
      typename abstract_interface<T>::type>::type result_type;

   typedef S const& first_argument_type;
   typedef T const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      result_type R(size(y) + 1);
      zero_all(R);
      set_element(R, 0, x);
      typename const_iterator<T>::type I = iterate(y);
      while (I)
      {
         set_element(R, I.index() + 1, *I);
         ++I;
      }
      return R;
   }
};

template <typename S, typename T, typename Sv, typename Tv, typename Ti>
struct DirectSumInterface<S, T, AnyScalar<Sv>, DENSE_VECTOR(Tv, Ti)>
{
   typedef typename make_value<typename Addition<Sv, Tv>::result_type>::type value_type;

   typedef typename make_vector_from_abstract<value_type,
      typename abstract_interface<T>::type>::type result_type;

   typedef S const& first_argument_type;
   typedef T const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      result_type R(size(y) + 1);
      set_element(R, 0, x);
      typename const_iterator<T>::type I = iterate(y);
      while (I)
      {
         set_element(R, I.index() + 1, *I);
         ++I;
      }
      return R;
   }
};

// specialization doesn't require Addition to be defined
template <typename S, typename T, typename Val, typename Ti>
struct DirectSumInterface<S, T, AnyScalar<Val>, LOCAL_VECTOR(Val, Ti)>
{
   typedef typename make_value<T>::type result_type;
   typedef S const& first_argument_type;
   typedef T const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      result_type R(size(y) + 1);
      zero_all(R);
      set_element(R, 0, x);
      typename const_iterator<T>::type I = iterate(y);
      while (I)
      {
         set_element(R, I.index() + 1, *I);
         ++I;
      }
      return R;
   }
};

template <typename S, typename T, typename Val, typename Ti>
struct DirectSumInterface<S, T, AnyScalar<Val>, DENSE_VECTOR(Val, Ti)>
{
   typedef typename make_value<T>::type result_type;
   typedef S const& first_argument_type;
   typedef T const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      result_type R(size(y) + 1);
      set_element(R, 0, x);
      typename const_iterator<T>::type I = iterate(y);
      while (I)
      {
         set_element(R, I.index() + 1, *I);
         ++I;
      }
      return R;
   }
};

} // namespace LinearAlgebra

#endif
