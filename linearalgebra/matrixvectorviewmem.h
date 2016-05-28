// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/matrixvectorviewmem.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER
//
// Defines the VectorView function for STRIDE_MATRIX variants.
//

#if !defined(MATRIXVECTORMEMVIEW_H_4328UVFDE8J3T897HBR78HPG98)
#define MATRIXVECTORMEMVIEW_H_4328UVFDE8J3T897HBR78HPG98

#include "matrixinterface.h"
#include "vectormemproxy.h"
#include "vector.h"

namespace LinearAlgebra
{

// for stride matrices, the best we can do is a vector of variable stride proxies.
// The loops here are very ugly; VectorMemProxy has reference semantics so we cannot
// use ordinary assignment, and there is no easy way to construct the vector of proxies
// correctly in place.

// const version

template <typename Vec, typename T, typename Ti>
struct VectorViewInterface<Vec,  Concepts::StrideRowMajorMatrix<T, Ti>>
{
   typedef Vec argument_type;
   typedef Vector<VectorMemProxy<T const, tagVariable> > result_type;
   result_type operator()(argument_type const& x) const
   {
      T const* d = data(x);
      int str1 = stride1(x);
      int str2 = stride2(x);
      size_t s1 = size1(x);
      size_t s2 = size2(x);
      result_type Result(s1);
      for (size_t i = 0; i < s1; ++i)
      {
         Result[i].~VectorMemProxy<T const, tagVariable>();
         new (&Result[i]) VectorMemProxy<T const, tagVariable>(d + i*str1, s2, str2);
      }
      return Result;
   }
};

template <typename Vec, typename T, typename Ti>
struct VectorViewInterface<Vec,  Concepts::StrideColMajorMatrix<T, Ti>>
{
   typedef Vec argument_type;
   typedef Vector<VectorMemProxy<T const, tagVariable> > result_type;
   result_type operator()(argument_type const& x) const
   {
      T const* d = data(x);
      int str1 = stride1(x);
      int str2 = stride2(x);
      size_t s1 = size1(x);
      size_t s2 = size2(x);
      result_type Result(s2);
      for (size_t i = 0; i < s2; ++i)
      {
         Result[i].~VectorMemProxy<T const, tagVariable>();
         new (&Result[i]) VectorMemProxy<T const, tagVariable>(d + i*str2, s1, str1);
      }
      return Result;
   }
};

// non-const version

template <typename Vec, typename T, typename Ti>
struct VectorViewInterface<Vec&, Concepts::StrideRowMajorMatrix<T, Ti>>
{
   typedef Vec argument_type;
   typedef Vector<VectorMemProxy<T, tagVariable> > result_type;
   result_type operator()(argument_type& x) const
   {
      T* d = data(x);
      int str1 = stride1(x);
      int str2 = stride2(x);
      size_t s1 = size1(x);
      size_t s2 = size2(x);
      result_type Result(s1);
      for (size_t i = 0; i < s1; ++i)
      {
         Result[i].~VectorMemProxy<T, tagVariable>();
         new (&Result[i]) VectorMemProxy<T, tagVariable>(d + i*str1, s2, str2);
      }
      return Result;
   }
};

template <typename Vec, typename T, typename Ti>
struct VectorViewInterface<Vec&, Concepts::StrideColMajorMatrix<T, Ti>>
{
   typedef Vec argument_type;
   typedef Vector<VectorMemProxy<T, tagVariable> > result_type;
   result_type operator()(argument_type& x) const
   {
      T* d = data(x);
      int str1 = stride1(x);
      int str2 = stride2(x);
      size_t s1 = size1(x);
      size_t s2 = size2(x);
      result_type Result(s2);
      for (size_t i = 0; i < s2; ++i)
      {
         Result[i].~VectorMemProxy<T, tagVariable>();
         new (&Result[i]) VectorMemProxy<T, tagVariable>(d + i*str2, s1, str1);
      }
      return Result;
   }
};

// for contiguous matrices, we can have a vector of stride 1 proxies
template <typename Vec, typename T, typename Ti>
struct VectorViewInterface<Vec,  Concepts::ContiguousRowMajorMatrix<T, Ti>>
{
   typedef Vec argument_type;
   typedef Vector<VectorMemProxy<T const> > result_type;
   result_type operator()(argument_type const& x) const
   {
      T const* d = data(x);
      int str1 = stride1(x);
      size_t s1 = size1(x);
      size_t s2 = size2(x);
      result_type Result(s1);
      for (size_t i = 0; i < s1; ++i)
      {
         Result[i].~VectorMemProxy<T const>();
         new (&Result[i]) VectorMemProxy<T const>(d + i*str1, s2);
      }
      return Result;
   }
};

template <typename Vec, typename T, typename Ti>
struct VectorViewInterface<Vec, Concepts::ContiguousColMajorMatrix<T, Ti>>
{
   typedef Vec argument_type;
   typedef Vector<VectorMemProxy<T const> > result_type;
   result_type operator()(argument_type const& x) const
   {
      T const* d = data(x);
      int str2 = stride2(x);
      size_t s1 = size1(x);
      size_t s2 = size2(x);
      result_type Result(s2);
      for (size_t i = 0; i < s2; ++i)
      {
         Result[i].~VectorMemProxy<T const>();
         new (&Result[i]) VectorMemProxy<T const>(d + i*str2, s1);
      }
      return Result;
   }
};

// non-const

template <typename Vec, typename T, typename Ti>
struct VectorViewInterface<Vec&,  Concepts::ContiguousRowMajorMatrix<T, Ti>>
{
   typedef Vec argument_type;
   typedef Vector<VectorMemProxy<T> > result_type;
   result_type operator()(argument_type& x) const
   {
      T* d = data(x);
      int str1 = stride1(x);
      size_t s1 = size1(x);
      size_t s2 = size2(x);
      result_type Result(s1);
      for (size_t i = 0; i < s1; ++i)
      {
         Result[i].~VectorMemProxy<T>();
         new (&Result[i]) VectorMemProxy<T>(d + i*str1, s2);
      }
      return Result;
   }
};

template <typename Vec, typename T, typename Ti>
struct VectorViewInterface<Vec&, Concepts::ContiguousColMajorMatrix<T, Ti>>
{
   typedef Vec argument_type;
   typedef Vector<VectorMemProxy<T> > result_type;
   result_type operator()(argument_type& x) const
   {
      T* d = data(x);
      int str2 = stride2(x);
      size_t s1 = size1(x);
      size_t s2 = size2(x);
      result_type Result(s2);
      for (size_t i = 0; i < s2; ++i)
      {
         Result[i].~VectorMemProxy<T>();
         new (&Result[i]) VectorMemProxy<T>(d + i*str2, s1);
      }
      return Result;
   }
};

} // namespace

#endif
