// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/matrixdirectproduct.h
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
  matrixdirectproduct.h

  direct_product function for matrices.

  Created 2005-03-02 Ian McCulloch

  FIXME:
  This is a very rough attempt - no expression templates, supports dense matrices only.

*/

#if !defined(MATRIXDIRECTPRODUCT_H_JCIFHUHUHLIURHGLO9H)
#define MATRIXDIRECTPRODUCT_H_JCIFHUHUHLIURHGLO9H

#include "matrixoperationsbase.h"
#include "index.h"
#include "matrix.h"
#include "matrixtransform.h"

namespace LinearAlgebra
{

template <typename T1, typename T2, typename Nested,
          typename T1v, typename T1o, typename T1i,
          typename T2v, typename T2o, typename T2i>
struct MatrixDirectProduct<T1, T2, Nested,
                           Concepts::DenseMatrix<T1v, T1o, T1i>,
                           Concepts::StrideMatrix<T2v, T2o, T2i>>
{
   typedef T1 const& first_argument_type;
   typedef T2 const& second_argument_type;
   typedef Nested const& third_argument_type;
   typedef typename make_value<typename Nested::result_type>::type scalar_type;
   typedef Matrix<scalar_type, T2o> result_type;

   result_type operator()(T1 const& m1, T2 const& m2, Nested const& f) const
   {
      size_type m1s1 = size1(m1);
      size_type m1s2 = size2(m1);
      size_type m2s1 = size1(m2);
      size_type m2s2 = size2(m2);
      result_type Res(m1s1 * m2s1, m1s2 * m2s2);
      for (size_type i = 0; i < m1s1; ++i)
      {
         for (size_type j = 0; j < m1s2; ++j)
         {
            //      TRACE("here")(tracer::typeid_name(*this));
            assign(project(Res, Range(i * m2s1, (i+1) * m2s1), Range(j * m2s2, (j+1) * m2s2)),
                   left_scalar_prod(get_element(m1, i,j), m2, f));
         }
      }
      return Res;
   }

   result_type operator()(T1 const& m1, T2 const& m2) const
   {
      return this->operator()(m1, m2, Nested());
   }
};

// FIXME: this should be sparse

template <typename T1, typename T2, typename Nested,
          typename T1v, typename T1o, typename T1i,
          typename T2v, typename T2i>
struct MatrixDirectProduct<T1, T2, Nested,
                           Concepts::DenseMatrix<T1v, T1o, T1i>,
                           Concepts::SparseMatrix<T2v, T2i>>
{
   typedef T1 const& first_argument_type;
   typedef T2 const& second_argument_type;
   typedef Nested const& third_argument_type;
   typedef typename make_value<typename Nested::result_type>::type scalar_type;
   typedef Matrix<scalar_type, T1o> result_type;

   result_type operator()(T1 const& m1, T2 const& m2, Nested const& f) const
   {
      size_type m1s1 = size1(m1);
      size_type m1s2 = size2(m1);
      size_type m2s1 = size1(m2);
      size_type m2s2 = size2(m2);
      result_type Res(m1s1 * m2s1, m1s2 * m2s2);
      for (size_type i = 0; i < m1s1; ++i)
      {
         for (size_type j = 0; j < m1s2; ++j)
         {
            //      TRACE("here")(tracer::typeid_name(*this));
            assign(project(Res, Range(i * m2s1, (i+1) * m2s1), Range(j * m2s2, (j+1) * m2s2)),
                   left_scalar_prod(get_element(m1, i,j), m2, f));
         }
      }
      return Res;
   }

   result_type operator()(T1 const& m1, T2 const& m2) const
   {
      return this->operator()(m1, m2, Nested());
   }
};

// FIXME: this should be sparse, or even block-sparse

template <typename T1, typename T2, typename Nested,
          typename T1v, typename T1i,
          typename T2v, typename T2o, typename T2i>
struct MatrixDirectProduct<T1, T2, Nested,
                           Concepts::SparseMatrix<T1v, T1i>,
                           Concepts::DenseMatrix<T2v, T2o, T2i>>
{
   typedef T1 const& first_argument_type;
   typedef T2 const& second_argument_type;
   typedef Nested const& third_argument_type;
   typedef typename make_value<typename Nested::result_type>::type scalar_type;
   typedef Matrix<scalar_type, T2o> result_type;

   result_type operator()(T1 const& m1, T2 const& m2, Nested const& f) const
   {
      size_type m1s1 = size1(m1);
      size_type m1s2 = size2(m1);
      size_type m2s1 = size1(m2);
      size_type m2s2 = size2(m2);
      result_type Res(m1s1 * m2s1, m1s2 * m2s2);
      zero_all(Res);
      typename const_iterator<T1>::type I = iterate(m1);
      while (I)
      {
         typename const_inner_iterator<T1>::type J = iterate(I);
         while (J)
         {
            size_type i = J.index1();
            size_type j = J.index2();

            assign(project(Res, Range(i * m2s1, (i+1) * m2s1), Range(j * m2s2, (j+1) * m2s2)),
                   left_scalar_prod(*J, m2, f));
            ++J;
         }
         ++I;
      }
      return Res;
   }

   result_type operator()(T1 const& m1, T2 const& m2) const
   {
      return this->operator()(m1, m2, Nested());
   }
};


template <typename T1, typename T2, typename Nested,
          typename T1v, typename T1o, typename T1i,
          typename T2v, typename T2o, typename T2i>
struct MatrixDirectProduct<T1, T2, Nested,
                           Concepts::DenseMatrix<T1v, T1o, T1i>,
                           Concepts::DenseMatrix<T2v, T2o, T2i>>
   : MatrixDirectProduct<T1, Matrix<T2v, T2o>, Nested> {};

template <typename T1, typename T2, typename Nested,
          typename T1i,
          typename T2v, typename T2i>
struct MatrixDirectProduct<T1, T2, Nested,
                           AnyScalar<T1i>,
                           Concepts::AnyMatrix<T2v, T2i>>
: ScalarMatrixMultiplication<T1, T2, Nested> {};

template <typename T1, typename T2, typename Nested,
          typename T1v, typename T1i,
          typename T2i>
struct MatrixDirectProduct<T1, T2, Nested,
                           Concepts::AnyMatrix<T1v, T1i>,
                           AnyScalar<T2i> >
: MatrixScalarMultiplication<T1, T2, Nested> {};

} // namespace LinearAlgebra

#endif
