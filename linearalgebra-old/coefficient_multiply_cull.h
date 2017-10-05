// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/coefficient_multiply_cull.h
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

#if !defined(MPTOOLKIT_LINEARALGEBRA_COEFFICIENT_MULTIPLY_CULL_H)
#define MPTOOLKIT_LINEARALGEBRA_COEFFICIENT_MULTIPLY_CULL_H

#include "coefficient_multiply.h"

namespace LinearAlgebra
{

//
// assign_coefficient_cull_product2
//

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSInterface = typename interface<LHS>::type,
          typename M1Interface = typename interface<M1>::type,
          typename M2Interface = typename interface<M2>::type>
struct AssignCoefficientCull_Product2 {};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float>
inline
void assign_coefficient_cull_product2(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
{
   AssignCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float>::apply(lhs, m1, m2, cf, f, Tol);
}

template <typename LHS, typename M1, typename M2,
          typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSi,
          typename M1v, typename M1i,
          typename M2v, typename M2i>
struct AssignCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                      Concepts::MatrixExpression<LHSv, LHSi>,
                      Concepts::MatrixExpression<M1v, M1i>,
                      Concepts::MatrixExpression<M2v, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      assign_coefficient_cull_product2(lhs, eval_expression(m1), eval_expression(m2), cf, f, Tol);
   }
};

// assign dense

template <typename LHS, typename M1, typename M2,
          typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSOrient, typename LHSi,
          typename M1v, typename M1Orient, typename M1i,
          typename M2v, typename M2Orient, typename M2i>
struct AssignCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
                      Concepts::DenseMatrix<M1v, M1Orient, M1i>,
                      Concepts::DenseMatrix<M2v, M2Orient, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      Matrix<M1v, RowMajor> Temp1(m1);
      Matrix<M2v, ColMajor> Temp2(m2);
      Matrix<LHSv, RowMajor> TempL(size1(m1), size2(m2));
      assign_coefficient_cull_product2(TempL, Temp1, Temp2, cf, f, Tol);
      assign(lhs, TempL);
   }
};

template <typename LHS, typename M1, typename M2,
          typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSOrient, typename LHSi,
          typename M1v, typename M1Orient, typename M1i,
          typename M2v, typename M2Orient, typename M2i>
struct AssignCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                      Concepts::StrideMatrix<LHSv, LHSOrient, LHSi>,
                      Concepts::DenseMatrix<M1v, M1Orient, M1i>,
                      Concepts::DenseMatrix<M2v, M2Orient, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      Matrix<M1v, RowMajor> Temp1(m1);
      Matrix<M2v, ColMajor> Temp2(m2);
      assign_coefficient_cull_product2(lhs, Temp1, Temp2, cf, f, Tol);
   }
};

// a partial specialization for LHS = RowMajor * ColMajor
// would be rude here - it would force specializations on the
// LHS or the value_type to also follow suit.
// Instead, we forward to AssignCoefficientCull_Product2_StrideStrideStride
// and specialize further there.

template <typename LHS, typename LHSOrient,
          typename M1, typename M1Orient,
          typename M2, typename M2Orient,
          typename CF,
          typename Nested, typename Float>
struct AssignCoefficientCull_Product2_StrideStrideStride;

template <typename LHS, typename LHSOrient,
          typename M1, typename M2,
          typename CF,
          typename Nested, typename Float>
struct AssignCoefficientCull_Product2_StrideStrideStride<LHS, LHSOrient, M1, RowMajor, M2, ColMajor,
                                                         CF, Nested, Float>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      typename const_iterator<M1>::type I = iterate(m1);
      while (I)
      {
         typename const_iterator<M2>::type J = iterate(m2);
         while (J)
         {
            set_element(lhs, I.index(), J.index(),
                        coefficient_parallel_prod_cull(*I, *J, bindij(cf, I.index(), J.index()), f, Tol));
            ++J;
         }
         ++I;
      }
   }
};

template <typename LHS, typename LHSOrient,
          typename M1, typename M2, typename CF, typename Nested, typename Float>
struct AssignCoefficientCull_Product2_StrideStrideStride<LHS, LHSOrient, M1, ColMajor, M2, ColMajor,
                                                         CF, Nested, Float>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      assign_coefficient_cull_product2(lhs, swap_sort_order(m1), m2, cf, f, Tol);
   }
};

template <typename LHS, typename LHSOrient,
          typename M1, typename M2, typename CF, typename Nested, typename Float>
struct AssignCoefficientCull_Product2_StrideStrideStride<LHS, LHSOrient, M1, RowMajor, M2, RowMajor,
                                                         CF, Nested, Float>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      assign_coefficient_cull_product2(lhs, m1, swap_sort_order(m2), cf, f, Tol);
   }
};

template <typename LHS, typename LHSOrient,
          typename M1, typename M2, typename CF, typename Nested, typename Float>
struct AssignCoefficientCull_Product2_StrideStrideStride<LHS, LHSOrient, M1, ColMajor, M2, RowMajor,
                                                         CF, Nested, Float>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested f, Float const& Tol)
   {
      assign_coefficient_cull_product2(lhs, swap_sort_order(m1), swap_sort_order(m2), cf, f, Tol);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSOrient, typename LHSi,
          typename M1v, typename M1Orient, typename M1i,
          typename M2v, typename M2Orient, typename M2i>
          struct AssignCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                                                Concepts::StrideMatrix<LHSv, LHSOrient, LHSi>,
                                                Concepts::StrideMatrix<M1v, M1Orient, M1i>,
                                                Concepts::StrideMatrix<M2v, M2Orient, M2i>>
          : AssignCoefficientCull_Product2_StrideStrideStride<LHS, LHSOrient,
                                                              M1, M1Orient,
                                                              M2, M2Orient,
                                                              CF, Nested, Float> {};

// assign sparse

template <typename LHS,
          typename M1, typename M1Orient,
          typename M2, typename M2Orient,
          typename CF, typename Nested, typename Float>
struct AssignCoefficientCull_Product2_Sparse;

template <typename LHS,
          typename M1, typename M2, typename CF, typename Nested, typename Float>
struct AssignCoefficientCull_Product2_Sparse<LHS, M1, RowMajor, M2, ColMajor, CF, Nested, Float>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      zero_all(lhs);
      typename const_iterator<M1>::type I = iterate(m1);
      while (I)
      {
         typename const_iterator<M2>::type J = iterate(m2);
         while (J)
         {
            add_element_cull(lhs, I.index(), J.index(),
                             coefficient_parallel_prod_cull(*I, *J,
                                                            bindij(cf, I.index(), J.index()), Tol, f), Tol);
            ++J;
         }
         ++I;
      }
   }
};

template <typename LHS,
          typename M1, typename M2, typename CF, typename Nested, typename Float>
struct AssignCoefficientCull_Product2_Sparse<LHS, M1, ColMajor, M2, ColMajor, CF, Nested, Float>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      assign_coefficient_cull_product2(lhs, Temp1, m2, cf, f, Tol);
   }
};

template <typename LHS,
          typename M1, typename M2, typename CF, typename Nested, typename Float>
struct AssignCoefficientCull_Product2_Sparse<LHS, M1, ColMajor, M2, RowMajor, CF, Nested, Float>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      assign_coefficient_cull_product2(lhs, Temp1, Temp2, cf, f, Tol);
   }
};

template <typename LHS,
          typename M1, typename M2, typename CF, typename Nested, typename Float>
struct AssignCoefficientCull_Product2_Sparse<LHS, M1, RowMajor, M2, RowMajor, CF, Nested, Float>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      assign_coefficient_cull_product2(lhs, m1, Temp2, cf, f, Tol);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSi,
          typename M1v, typename M1Orient, typename M1i,
          typename M2v, typename M2Orient, typename M2i>
struct AssignCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                      Concepts::LocalMatrix<LHSv, LHSi>,
                      Concepts::CompressedOuterMatrix<M1v, M1Orient, M1i>,
                      Concepts::CompressedOuterMatrix<M2v, M2Orient, M2i>>
: AssignCoefficientCull_Product2_Sparse<LHS, M1, M1Orient, M2, M2Orient, CF, Nested, Float> {};

// assign mixed sparse/dense

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSOrient, typename LHSi,
          typename M1v, typename M1i,
          typename M2v, typename M2i>
struct AssignCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
                      Concepts::CompressedOuterMatrix<M1v, RowMajor, M1i>,
                      Concepts::DenseMatrix<M2v, ColMajor, M2i>>
: AssignCoefficientCull_Product2_Sparse<LHS, M1, RowMajor, M2, ColMajor, CF, Nested, Float> {};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSOrient, typename LHSi,
          typename M1v, typename M1i,
          typename M2v, typename M2i>
struct AssignCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
                      Concepts::CompressedOuterMatrix<M1v, ColMajor, M1i>,
                      Concepts::DenseMatrix<M2v, ColMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      assign_coefficient_cull_product2(lhs, Temp1, m2, cf, f, Tol);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSOrient, typename LHSi,
          typename M1v, typename M1i,
          typename M2v, typename M2i>
struct AssignCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
                      Concepts::CompressedOuterMatrix<M1v, ColMajor, M1i>,
                      Concepts::DenseMatrix<M2v, RowMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      assign_coefficient_cull_product2(lhs, Temp1, swap_sort_order(m2), cf, f, Tol);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSOrient, typename LHSi,
          typename M1v, typename M1i,
          typename M2v, typename M2i>
struct AssignCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                                      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
                                      Concepts::CompressedOuterMatrix<M1v, RowMajor, M1i>,
                                      Concepts::DenseMatrix<M2v, RowMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      assign_coefficient_cull_product2(lhs, m1, swap_sort_order(m2), cf, f, Tol);
   }
};


template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSOrient, typename LHSi,
          typename M1v, typename M1i,
          typename M2v, typename M2i>
struct AssignCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                                      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
                                      Concepts::DenseMatrix<M1v, RowMajor, M1i>,
                                      Concepts::CompressedOuterMatrix<M2v, ColMajor, M2i>>
: AssignCoefficientCull_Product2_Sparse<LHS, M1, RowMajor, M2, ColMajor, CF, Nested, Float> {};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSOrient, typename LHSi,
          typename M1v, typename M1i,
          typename M2v, typename M2i>
struct AssignCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
                      Concepts::DenseMatrix<M1v, ColMajor, M1i>,
                      Concepts::CompressedOuterMatrix<M2v, ColMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      assign_coefficient_cull_product2(lhs, swap_sort_order(m1), m2, cf, f, Tol);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSOrient, typename LHSi,
          typename M1v, typename M1i,
          typename M2v, typename M2i>
 struct AssignCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
                      Concepts::DenseMatrix<M1v, ColMajor, M1i>,
                      Concepts::CompressedOuterMatrix<M2v, RowMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      assign_coefficient_cull_product2(lhs, swap_sort_order(m1), Temp2, cf, f, Tol);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSOrient, typename LHSi,
          typename M1v, typename M1i,
          typename M2v, typename M2i>
 struct AssignCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
                      Concepts::DenseMatrix<M1v, RowMajor, M1i>,
                      Concepts::CompressedOuterMatrix<M2v, RowMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      assign_coefficient_cull_product2(lhs, m1, Temp2, cf, f, Tol);
   }
};

// assign mixed sparse/diagonal

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSi,
          typename M1v, typename M1i,
          typename M2v, typename M2i>
struct AssignCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                                  Concepts::SparseMatrix<LHSv, LHSi>,
                                  Concepts::SparseMatrix<M1v, M1i>,
                                  Concepts::DiagonalMatrix<M2v, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      zero_all(lhs);
      typename const_iterator<M1>::type I = iterate(m1);
      while (I)
      {
         typename const_inner_iterator<M1>::type J = iterate(I);
         while (J)
         {
            add_element_cull(lhs, J.index1(), J.index2(),
                             cf(J.index1(), J.index2(), J.index2()) * f(*J, m2.diagonal()[J.index2()]), Tol);
            ++J;
         }
         ++I;
      }
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSi,
          typename M1v, typename M1i,
          typename M2v, typename M2i>
struct AssignCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                                  Concepts::SparseMatrix<LHSv, LHSi>,
                                  Concepts::DiagonalMatrix<M1v, M1i>,
                                  Concepts::SparseMatrix<M2v, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      zero_all(lhs);
      typename const_iterator<M2>::type I = iterate(m2);
      while (I)
      {
         typename const_inner_iterator<M2>::type J = iterate(I);
         while (J)
         {
            add_element_cull(lhs, J.index1(), J.index2(),
                             cf(J.index1(), J.index1(), J.index2()) * f(m1.diagonal()[J.index1()], *J), Tol);
            ++J;
         }
         ++I;
      }
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSi,
          typename M1v, typename M1i,
          typename M2v, typename M2i>
struct AssignCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                                  Concepts::SparseMatrix<LHSv, LHSi>,
                                  Concepts::DiagonalMatrix<M1v, M1i>,
                                  Concepts::DiagonalMatrix<M2v, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      zero_all(lhs);
      auto I = iterate_diagonal(m1);
      auto J = iterate_diagonal(m2);
      while (I)
      {
         add_element_cull(lhs, I.index(), I.index(),
                          cf(I.index(), I.index(), I.index()) * f(*I, *J), Tol);
         ++I;
         ++J;
      }
      DEBUG_CHECK(!J);
   }
};

//
// add_coefficient_cull_product2
//

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSInterface = typename interface<LHS>::type,
          typename M1Interface = typename interface<M1>::type,
          typename M2Interface = typename interface<M2>::type>
struct AddCoefficientCull_Product2 {};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float>
inline
void add_coefficient_cull_product2(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
{
   AddCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float>::apply(lhs, m1, m2, cf, f, Tol);
}

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSi,
          typename M1v, typename M1i,
          typename M2v, typename M2i>
struct AddCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                      Concepts::MatrixExpression<LHSv, LHSi>,
                      Concepts::MatrixExpression<M1v, M1i>,
                      Concepts::MatrixExpression<M2v, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      add_coefficient_cull_product2(lhs, eval_expression(m1), eval_expression(m2), cf, f, Tol);
   }
};

// add dense

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSOrient, typename LHSi,
          typename M1v, typename M1Orient, typename M1i,
          typename M2v, typename M2Orient, typename M2i>
struct AddCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
                      Concepts::DenseMatrix<M1v, M1Orient, M1i>,
                      Concepts::DenseMatrix<M2v, M2Orient, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      Matrix<M1v, RowMajor> Temp1(m1);
      Matrix<M2v, ColMajor> Temp2(m2);
      Matrix<LHSv, RowMajor> TempL(size1(m1), size2(m2));
      assign_coefficient_cull_product2(TempL, Temp1, Temp2, cf, f, Tol);
      add(lhs, TempL);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSOrient, typename LHSi,
          typename M1v, typename M1Orient, typename M1i,
          typename M2v, typename M2Orient, typename M2i>
struct AddCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                      Concepts::StrideMatrix<LHSv, LHSOrient, LHSi>,
                      Concepts::DenseMatrix<M1v, M1Orient, M1i>,
                      Concepts::DenseMatrix<M2v, M2Orient, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      Matrix<M1v, RowMajor> Temp1(m1);
      Matrix<M2v, ColMajor> Temp2(m2);
      add_coefficient_cull_product2(lhs, Temp1, Temp2, cf, f, Tol);
   }
};

// a partial specialization for LHS = RowMajor * ColMajor
// would be rude here - it would force specializations on the
// LHS or the value_type to also follow suit.
// Instead, we forward to AddCoefficientCull_Product2_StrideStrideStride
// and specialize further there.

template <typename LHS, typename LHSOrient,
          typename M1, typename M1Orient,
          typename M2, typename M2Orient,
          typename CF, typename Nested, typename Float>
struct AddCoefficientCull_Product2_StrideStrideStride;

template <typename LHS, typename LHSOrient,
          typename M1, typename M2,
          typename CF, typename Nested, typename Float>
struct AddCoefficientCull_Product2_StrideStrideStride<LHS, LHSOrient, M1, RowMajor, M2, ColMajor,
                                                      CF, Nested, Float>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      typename const_iterator<M1>::type I = iterate(m1);
      while (I)
      {
         typename const_iterator<M2>::type J = iterate(m2);
         while (J)
         {
            add_element_cull(lhs, I.index(), J.index(),
                             coefficient_parallel_prod_cull(*I, *J, bindij(cf, I.index(), J.index()), f, Tol), Tol);
            ++J;
         }
         ++I;
      }
   }
};

template <typename LHS, typename LHSOrient,
          typename M1, typename M2, typename CF, typename Nested, typename Float>
struct AddCoefficientCull_Product2_StrideStrideStride<LHS, LHSOrient, M1, ColMajor, M2, ColMajor,
                                                      CF, Nested, Float>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      add_coefficient_cull_product2(lhs, swap_sort_order(m1), m2, cf, f, Tol);
   }
};

template <typename LHS, typename LHSOrient,
          typename M1, typename M2, typename CF, typename Nested, typename Float>
struct AddCoefficientCull_Product2_StrideStrideStride<LHS, LHSOrient, M1, RowMajor, M2, RowMajor,
                                                      CF, Nested, Float>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      add_coefficient_cull_product2(lhs, m1, swap_sort_order(m2), cf, f, Tol);
   }
};

template <typename LHS, typename LHSOrient,
          typename M1, typename M2, typename CF, typename Nested, typename Float>
struct AddCoefficientCull_Product2_StrideStrideStride<LHS, LHSOrient, M1, ColMajor, M2, RowMajor,
                                                      CF, Nested, Float>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      add_coefficient_cull_product2(lhs, swap_sort_order(m1), swap_sort_order(m2), cf, f, Tol);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSOrient, typename LHSi,
          typename M1v, typename M1Orient, typename M1i,
          typename M2v, typename M2Orient, typename M2i>
struct AddCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                      Concepts::StrideMatrix<LHSv, LHSOrient, LHSi>,
                      Concepts::StrideMatrix<M1v, M1Orient, M1i>,
                      Concepts::StrideMatrix<M2v, M2Orient, M2i>>
: AddCoefficientCull_Product2_StrideStrideStride<LHS, LHSOrient, M1, M1Orient, M2, M2Orient,
                                                 CF, Nested, Float> {};


// add sparse

template <typename LHS,
          typename M1, typename M1Orient,
          typename M2, typename M2Orient, typename CF, typename Nested, typename Float>
struct AddCoefficientCull_Product2_Sparse;

template <typename LHS,
          typename M1, typename M2, typename CF, typename Nested, typename Float>
struct AddCoefficientCull_Product2_Sparse<LHS, M1, RowMajor, M2, ColMajor, CF, Nested, Float>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      typename const_iterator<M1>::type I = iterate(m1);
      while (I)
      {
         typename const_iterator<M2>::type J = iterate(m2);
         while (J)
         {
            add_element_check_if_zero(lhs, I.index(), J.index(),
                                      coefficient_parallel_prod_cull(*I, *J, bindij(cf, I.index(), J.index()), f, Tol));
            ++J;
         }
         ++I;
      }
   }
};

template <typename LHS,
          typename M1, typename M2, typename CF, typename Nested, typename Float>
struct AddCoefficientCull_Product2_Sparse<LHS, M1, ColMajor, M2, ColMajor, CF, Nested, Float>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      add_coefficient_cull_product2(lhs, Temp1, m2, cf, f, Tol);
   }
};

template <typename LHS,
          typename M1, typename M2, typename CF, typename Nested, typename Float>
struct AddCoefficientCull_Product2_Sparse<LHS, M1, ColMajor, M2, RowMajor, CF, Nested, Float>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      add_coefficient_cull_product2(lhs, Temp1, Temp2, cf, f, Tol);
   }
};

template <typename LHS,
          typename M1, typename M2, typename CF, typename Nested, typename Float>
struct AddCoefficientCull_Product2_Sparse<LHS, M1, RowMajor, M2, RowMajor, CF, Nested, Float>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      add_coefficient_cull_product2(lhs, m1, Temp2, cf, f, Tol);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSi,
          typename M1v, typename M1Orient, typename M1i,
          typename M2v, typename M2Orient, typename M2i>
struct AddCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                      Concepts::LocalMatrix<LHSv, LHSi>,
                      Concepts::CompressedOuterMatrix<M1v, M1Orient, M1i>,
                      Concepts::CompressedOuterMatrix<M2v, M2Orient, M2i>>
: AddCoefficientCull_Product2_Sparse<LHS, M1, M1Orient, M2, M2Orient, CF, Nested, Float> {};

// add mixed sparse/dense

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSOrient, typename LHSi,
          typename M1v, typename M1i,
          typename M2v, typename M2i>
struct AddCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
                      Concepts::CompressedOuterMatrix<M1v, RowMajor, M1i>,
                      Concepts::DenseMatrix<M2v, ColMajor, M2i>>
: AddCoefficientCull_Product2_Sparse<LHS, M1, RowMajor, M2, ColMajor, CF, Nested, Float> {};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSOrient, typename LHSi,
          typename M1v, typename M1i,
          typename M2v, typename M2i>
struct AddCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
                      Concepts::CompressedOuterMatrix<M1v, ColMajor, M1i>,
                      Concepts::DenseMatrix<M2v, ColMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      add_coefficient_cull_product2(lhs, Temp1, m2, cf, f, Tol);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSOrient, typename LHSi,
          typename M1v, typename M1i,
          typename M2v, typename M2i>
struct AddCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
                      Concepts::CompressedOuterMatrix<M1v, ColMajor, M1i>,
                      Concepts::DenseMatrix<M2v, RowMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      add_coefficient_cull_product2(lhs, Temp1, swap_sort_order(m2), cf, f, Tol);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSOrient, typename LHSi,
          typename M1v, typename M1i,
          typename M2v, typename M2i>
struct AddCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
                      Concepts::CompressedOuterMatrix<M1v, RowMajor, M1i>,
                      Concepts::DenseMatrix<M2v, RowMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      add_coefficient_cull_product2(lhs, m1, swap_sort_order(m2), cf, f, Tol);
   }
};


template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSOrient, typename LHSi,
          typename M1v, typename M1i,
          typename M2v, typename M2i>
struct AddCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
                      Concepts::DenseMatrix<M1v, RowMajor, M1i>,
                      Concepts::CompressedOuterMatrix<M2v, ColMajor, M2i>>
: AddCoefficientCull_Product2_Sparse<LHS, M1, RowMajor, M2, ColMajor, CF, Nested, Float> {};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSOrient, typename LHSi,
          typename M1v, typename M1i,
          typename M2v, typename M2i>
struct AddCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
                      Concepts::DenseMatrix<M1v, ColMajor, M1i>,
                      Concepts::CompressedOuterMatrix<M2v, ColMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      add_coefficient_cull_product2(lhs, swap_sort_order(m1), m2, cf, f, Tol);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSOrient, typename LHSi,
          typename M1v, typename M1i,
          typename M2v, typename M2i>
struct AddCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
                      Concepts::DenseMatrix<M1v, ColMajor, M1i>,
                      Concepts::CompressedOuterMatrix<M2v, RowMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      add_coefficient_cull_product2(lhs, swap_sort_order(m1), Temp2, cf, f, Tol);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSOrient, typename LHSi,
          typename M1v, typename M1i,
          typename M2v, typename M2i>
struct AddCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
                      Concepts::DenseMatrix<M1v, RowMajor, M1i>,
                      Concepts::CompressedOuterMatrix<M2v, RowMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      add_coefficient_cull_product2(lhs, m1, Temp2, cf, f, Tol);
   }
};

// add mixed sparse/diagonal

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSi,
          typename M1v, typename M1i,
          typename M2v, typename M2i>
struct AddCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                               Concepts::SparseMatrix<LHSv, LHSi>,
                               Concepts::SparseMatrix<M1v, M1i>,
                               Concepts::DiagonalMatrix<M2v, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      typename const_iterator<M1>::type I = iterate(m1);
      while (I)
      {
         typename const_inner_iterator<M1>::type J = iterate(I);
         while (J)
         {
            add_element_cull(lhs, J.index1(), J.index2(),
                             cf(J.index1(), J.index2(), J.index2()) * f(*J, m2.diagonal()[J.index2()]), Tol);
            ++J;
         }
         ++I;
      }
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSi,
          typename M1v, typename M1i,
          typename M2v, typename M2i>
struct AddCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                               Concepts::SparseMatrix<LHSv, LHSi>,
                               Concepts::DiagonalMatrix<M1v, M1i>,
                               Concepts::SparseMatrix<M2v, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      typename const_iterator<M2>::type I = iterate(m2);
      while (I)
      {
         typename const_inner_iterator<M2>::type J = iterate(I);
         while (J)
         {
            add_element_cull(lhs, J.index1(), J.index2(),
                             cf(J.index1(), J.index1(), J.index2()) * f(m1.diagonal()[J.index1()], *J), Tol);
            ++J;
         }
         ++I;
      }
   }
};

//
// subtract_coefficient_cull_product2
//

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSInterface = typename interface<LHS>::type,
          typename M1Interface = typename interface<M1>::type,
          typename M2Interface = typename interface<M2>::type>
struct SubtractCoefficientCull_Product2 {};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float>
inline
void subtract_coefficient_cull_product2(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
{
   SubtractCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float>::apply(lhs, m1, m2, cf, f, Tol);
}

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSi,
          typename M1v, typename M1i,
          typename M2v, typename M2i>
struct SubtractCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                      Concepts::MatrixExpression<LHSv, LHSi>,
                      Concepts::MatrixExpression<M1v, M1i>,
                      Concepts::MatrixExpression<M2v, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      subtract_coefficient_cull_product2(lhs, eval_expression(m1), eval_expression(m2), cf, f, Tol);
   }
};

// subtract dense

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSOrient, typename LHSi,
          typename M1v, typename M1Orient, typename M1i,
          typename M2v, typename M2Orient, typename M2i>
struct SubtractCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
                      Concepts::DenseMatrix<M1v, M1Orient, M1i>,
                      Concepts::DenseMatrix<M2v, M2Orient, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      Matrix<M1v, RowMajor> Temp1(m1);
      Matrix<M2v, ColMajor> Temp2(m2);
      Matrix<LHSv, RowMajor> TempL(size1(m1), size2(m2));
      assign_coefficient_cull_product2(TempL, Temp1, Temp2, cf, f, Tol);
      subtract(lhs, TempL);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSOrient, typename LHSi,
          typename M1v, typename M1Orient, typename M1i,
          typename M2v, typename M2Orient, typename M2i>
struct SubtractCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                      Concepts::StrideMatrix<LHSv, LHSOrient, LHSi>,
                      Concepts::DenseMatrix<M1v, M1Orient, M1i>,
                      Concepts::DenseMatrix<M2v, M2Orient, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      Matrix<M1v, RowMajor> Temp1(m1);
      Matrix<M2v, ColMajor> Temp2(m2);
      subtract_coefficient_cull_product2(lhs, Temp1, Temp2, cf, f, Tol);
   }
};

// a partial specialization for LHS = RowMajor * ColMajor
// would be rude here - it would force specializations on the
// LHS or the value_type to also follow suit.
// Instead, we forward to SubtractCoefficientCull_Product2_StrideStrideStride
// and specialize further there.

template <typename LHS, typename LHSOrient,
          typename M1, typename M1Orient,
          typename M2, typename M2Orient,
          typename CF, typename Nested, typename Float>
struct SubtractCoefficientCull_Product2_StrideStrideStride;

template <typename LHS, typename LHSOrient,
          typename M1, typename M2,
          typename CF, typename Nested, typename Float>
struct SubtractCoefficientCull_Product2_StrideStrideStride<LHS, LHSOrient, M1, RowMajor, M2, ColMajor,
                                                           CF, Nested, Float>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      typename const_iterator<M1>::type I = iterate(m1);
      while (I)
      {
         typename const_iterator<M2>::type J = iterate(m2);
         while (J)
         {
            subtract_element(lhs, I.index(), J.index(),
                             coefficient_parallel_prod_cull(*I, *J, bindij(cf, I.index(), J.index()), f, Tol));
            ++J;
         }
         ++I;
      }
   }
};

template <typename LHS, typename LHSOrient,
          typename M1, typename M2, typename CF, typename Nested, typename Float>
struct SubtractCoefficientCull_Product2_StrideStrideStride<LHS, LHSOrient, M1, ColMajor, M2, ColMajor,
                                                           CF, Nested, Float>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      subtract_coefficient_cull_product2(lhs, swap_sort_order(m1), m2, cf, f, Tol);
   }
};

template <typename LHS, typename LHSOrient,
          typename M1, typename M2, typename CF, typename Nested, typename Float>
struct SubtractCoefficientCull_Product2_StrideStrideStride<LHS, LHSOrient, M1, RowMajor, M2, RowMajor,
                                                           CF, Nested, Float>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      subtract_coefficient_cull_product2(lhs, m1, swap_sort_order(m2), cf, f, Tol);
   }
};

template <typename LHS, typename LHSOrient,
          typename M1, typename M2, typename CF, typename Nested, typename Float>
struct SubtractCoefficientCull_Product2_StrideStrideStride<LHS, LHSOrient, M1, ColMajor, M2, RowMajor,
                                                           CF, Nested, Float>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      subtract_coefficient_cull_product2(lhs, swap_sort_order(m1), swap_sort_order(m2), cf, f, Tol);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSOrient, typename LHSi,
          typename M1v, typename M1Orient, typename M1i,
          typename M2v, typename M2Orient, typename M2i>
struct SubtractCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                      Concepts::StrideMatrix<LHSv, LHSOrient, LHSi>,
                      Concepts::StrideMatrix<M1v, M1Orient, M1i>,
                      Concepts::StrideMatrix<M2v, M2Orient, M2i>>
: SubtractCoefficientCull_Product2_StrideStrideStride<LHS, LHSOrient, M1, M1Orient, M2, M2Orient,
                                                      CF, Nested, Float> {};

// subtract sparse

template <typename LHS,
          typename M1, typename M1Orient,
          typename M2, typename M2Orient, typename CF, typename Nested, typename Float>
struct SubtractCoefficientCull_Product2_Sparse;

template <typename LHS,
          typename M1, typename M2, typename CF, typename Nested, typename Float>
struct SubtractCoefficientCull_Product2_Sparse<LHS, M1, RowMajor, M2, ColMajor, CF, Nested, Float>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      typename const_iterator<M1>::type I = iterate(m1);
      while (I)
      {
         typename const_iterator<M2>::type J = iterate(m2);
         while (J)
         {
            subtract_element_check_if_zero(lhs, I.index(), J.index(),
                                           coefficient_parallel_prod_cull(*I, *J,
                                                                          bindij(cf, I.index(), J.index()),
                                                                          f, Tol));
            ++J;
         }
         ++I;
      }
   }
};

template <typename LHS,
          typename M1, typename M2, typename CF, typename Nested, typename Float>
struct SubtractCoefficientCull_Product2_Sparse<LHS, M1, ColMajor, M2, ColMajor, CF, Nested, Float>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      subtract_coefficient_cull_product2(lhs, Temp1, m2, cf, f, Tol);
   }
};

template <typename LHS,
          typename M1, typename M2, typename CF, typename Nested, typename Float>
struct SubtractCoefficientCull_Product2_Sparse<LHS, M1, ColMajor, M2, RowMajor, CF, Nested, Float>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      subtract_coefficient_cull_product2(lhs, Temp1, Temp2, cf, f, Tol);
   }
};

template <typename LHS,
          typename M1, typename M2, typename CF, typename Nested, typename Float>
struct SubtractCoefficientCull_Product2_Sparse<LHS, M1, RowMajor, M2, RowMajor, CF, Nested, Float>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      subtract_coefficient_cull_product2(lhs, m1, Temp2, cf, f, Tol);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSi,
          typename M1v, typename M1Orient, typename M1i,
          typename M2v, typename M2Orient, typename M2i>
struct SubtractCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                      Concepts::LocalMatrix<LHSv, LHSi>,
                      Concepts::CompressedOuterMatrix<M1v, M1Orient, M1i>,
                      Concepts::CompressedOuterMatrix<M2v, M2Orient, M2i>>
: SubtractCoefficientCull_Product2_Sparse<LHS, M1, M1Orient, M2, M2Orient, CF, Nested, Float> {};


// subtract mixed sparse/dense

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSOrient, typename LHSi,
          typename M1v, typename M1i,
          typename M2v, typename M2i>
struct SubtractCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
                      Concepts::CompressedOuterMatrix<M1v, RowMajor, M1i>,
                      Concepts::DenseMatrix<M2v, ColMajor, M2i>>
: SubtractCoefficientCull_Product2_Sparse<LHS, M1, RowMajor, M2, ColMajor, CF, Nested, Float> {};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSOrient, typename LHSi,
          typename M1v, typename M1i,
          typename M2v, typename M2i>
struct SubtractCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
                      Concepts::CompressedOuterMatrix<M1v, ColMajor, M1i>,
                      Concepts::DenseMatrix<M2v, ColMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      subtract_coefficient_cull_product2(lhs, Temp1, m2, cf, f, Tol);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSOrient, typename LHSi,
          typename M1v, typename M1i,
          typename M2v, typename M2i>
struct SubtractCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
                      Concepts::CompressedOuterMatrix<M1v, ColMajor, M1i>,
                      Concepts::DenseMatrix<M2v, RowMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      subtract_coefficient_cull_product2(lhs, Temp1, swap_sort_order(m2), cf, f, Tol);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSOrient, typename LHSi,
          typename M1v, typename M1i,
          typename M2v, typename M2i>
struct SubtractCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                                    Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
                                    Concepts::CompressedOuterMatrix<M1v, RowMajor, M1i>,
                                    Concepts::DenseMatrix<M2v, RowMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      subtract_coefficient_cull_product2(lhs, m1, swap_sort_order(m2), cf, f, Tol);
   }
};


template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSOrient, typename LHSi,
          typename M1v, typename M1i,
          typename M2v, typename M2i>
struct SubtractCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                                        Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
                                        Concepts::DenseMatrix<M1v, RowMajor, M1i>,
                                        Concepts::CompressedOuterMatrix<M2v, ColMajor, M2i>>
: SubtractCoefficientCull_Product2_Sparse<LHS, M1, RowMajor, M2, ColMajor, CF, Nested, Float> {};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSOrient, typename LHSi,
          typename M1v, typename M1i,
          typename M2v, typename M2i>
struct SubtractCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                                        Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
                                        Concepts::DenseMatrix<M1v, ColMajor, M1i>,
                                        Concepts::CompressedOuterMatrix<M2v, ColMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      subtract_coefficient_cull_product2(lhs, swap_sort_order(m1), m2, cf, f, Tol);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSOrient, typename LHSi,
          typename M1v, typename M1i,
          typename M2v, typename M2i>
struct SubtractCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                                        Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
                                        Concepts::DenseMatrix<M1v, ColMajor, M1i>,
                                        Concepts::CompressedOuterMatrix<M2v, RowMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      subtract_coefficient_cull_product2(lhs, swap_sort_order(m1), Temp2, cf, f, Tol);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSOrient, typename LHSi,
          typename M1v, typename M1i,
          typename M2v, typename M2i>
struct SubtractCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
                      Concepts::DenseMatrix<M1v, RowMajor, M1i>,
                      Concepts::CompressedOuterMatrix<M2v, RowMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      subtract_coefficient_cull_product2(lhs, m1, Temp2, cf, f, Tol);
   }
};

// subtract mixed sparse/diagonal

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSi,
          typename M1v, typename M1i,
          typename M2v, typename M2i>
struct SubtractCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                               Concepts::SparseMatrix<LHSv, LHSi>,
                               Concepts::SparseMatrix<M1v, M1i>,
                               Concepts::DiagonalMatrix<M2v, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      typename const_iterator<M1>::type I = iterate(m1);
      while (I)
      {
         typename const_inner_iterator<M1>::type J = iterate(I);
         while (J)
         {
            subtract_element_cull(lhs, J.index1(), J.index2(),
                                  cf(J.index1(), J.index2(), J.index2()) * f(*J, m2.diagonal()[J.index2()]), Tol);
            ++J;
         }
         ++I;
      }
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested, typename Float,
          typename LHSv, typename LHSi,
          typename M1v, typename M1i,
          typename M2v, typename M2i>
struct SubtractCoefficientCull_Product2<LHS, M1, M2, CF, Nested, Float,
                                                  Concepts::SparseMatrix<LHSv, LHSi>,
                                                  Concepts::DiagonalMatrix<M1v, M1i>,
                                                  Concepts::SparseMatrix<M2v, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f, Float const& Tol)
   {
      typename const_iterator<M2>::type I = iterate(m2);
      while (I)
      {
         typename const_inner_iterator<M2>::type J = iterate(I);
         while (J)
         {
            subtract_element_cull(lhs, J.index1(), J.index2(),
                                  cf(J.index1(), J.index1(), J.index2()) * f(m1.diagonal()[J.index1()], *J), Tol);
            ++J;
         }
         ++I;
      }
   }
};

} // namespace LinearAlgebra

#endif
