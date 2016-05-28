// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/matrixproductoperations.h
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

/*
   Created 2004-05-10 by Ian McCulloch

   The 'primitive' functions for matrix-matrix multiply.

   assign_product2(lhs, m1, m2)               implements lhs = prod(m1, m2)
   assign_scaled_product2(lhs, a, m1, m2)
      implements lhs = a * prod(m1, m2)  [or prod(a * m1, m2) etc]

   assign_product3(lhs, m1, m2, m3)     
      implements lhs = prod(m1, prod(m2, m3)) = prod(prod(m1, m2), m3)
      for matrices that are associative

   assign_scaled_product3(lhs, a, m1, m2, m3) 
      implements lhs = a * prod(m1, prod(m2, m3)) [and varations]

   add_product2(lhs, m1, m2)                  implements lhs += prod(m1, m2)
   add_scaled_product2(lhs, a, m1, m2)        implements lhs += a * prod(m1, m2)
   add_product3(lhs, m1, m2, m3)              implements lhs += prod(m1, prod(m2, m3))
   add_scaled_product3(lhs, a, m1, m2, m3)    implements lhs += a * prod(m1, prod(m2, m3))

   sub_product2(lhs, m1, m2)                  implements lhs -= prod(m1, m2)
   sub_scaled_product2(lhs, a, m1, m2)        implements lhs -= a * prod(m1, m2)
   sub_product3(lhs, m1, m2, m3)              implements lhs -= prod(m1, prod(m2, m3))
   sub_scaled_product3(lhs, a, m1, m2, m3)    implements lhs -= a * prod(m1, prod(m2, m3))
*/

#if !defined(MATRIXPRODUCT_H_DJSHC84YT789YUFP89RHJT89PJHFP8943)
#define MATRIXPRODUCT_H_DJSHC84YT789YUFP89RHJT89PJHFP8943

#include "matrixfwd.h"

namespace LinearAlgebra
{

//
// assign_product2
//

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSInterface = typename interface<LHS>::type,
	  typename M1Interface = typename interface<M1>::type,
	  typename M2Interface = typename interface<M2>::type>
struct AssignProduct2 {};

template <typename LHS, typename M1, typename M2, typename Nested>
inline
void assign_product2(LHS& lhs, M1 const& m1, M2 const& m2, Nested f)
{
   AssignProduct2<LHS, M1, M2, Nested>::apply(lhs, m1, m2, f);
}

template <typename LHS, typename M1, typename M2, typename Nested>
inline
typename boost::enable_if<is_mutable_proxy<LHS>, void>::type
assign_product2(LHS const& lhs, M1 const& m1, M2 const& m2, Nested f)
{
   AssignProduct2<LHS, M1, M2, Nested>::apply(const_cast<LHS&>(lhs), m1, m2, f);
}



template <typename LHS, typename M1, typename M2>
inline
void assign_product2(LHS& lhs, M1 const& m1, M2 const& m2)
{
   AssignProduct2<LHS, M1, M2, Multiplication<typename interface<M1>::value_type, 
      typename interface<M2>::value_type> >::
      apply(lhs, m1, m2, 
            Multiplication<typename interface<M1>::value_type,
                           typename interface<M2>::value_type>());
}

template <typename LHS, typename M1, typename M2,
          typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AssignProduct2<LHS, M1, M2, Nested,
		      Concepts::MatrixExpression<LHSv, LHSi>,
		      Concepts::MatrixExpression<M1v, M1i>,
		      Concepts::MatrixExpression<M2v, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested f)
   {
      assign_product2(lhs, eval_expression(m1), eval_expression(m2), f);
   }
};

// assign dense

template <typename LHS, typename M1, typename M2,
          typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1Orient, typename M1i,
	  typename M2v, typename M2Orient, typename M2i>
struct AssignProduct2<LHS, M1, M2, Nested,
		      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
		      Concepts::DenseMatrix<M1v, M1Orient, M1i>,
		      Concepts::DenseMatrix<M2v, M2Orient, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested f)
   {
      Matrix<M1v, RowMajor> Temp1(m1);
      Matrix<M2v, ColMajor> Temp2(m2);
      Matrix<LHSv, RowMajor> TempL(size1(m1), size2(m2));
      assign_product2(TempL, Temp1, Temp2, f);
      assign(lhs, TempL);
   }
};

template <typename LHS, typename M1, typename M2,
          typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1Orient, typename M1i,
	  typename M2v, typename M2Orient, typename M2i>
struct AssignProduct2<LHS, M1, M2, Nested,
		      Concepts::StrideMatrix<LHSv, LHSOrient, LHSi>,
		      Concepts::DenseMatrix<M1v, M1Orient, M1i>,
		      Concepts::DenseMatrix<M2v, M2Orient, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested f)
   {
      Matrix<M1v, RowMajor> Temp1(m1);
      Matrix<M2v, ColMajor> Temp2(m2);
      assign_product2(lhs, Temp1, Temp2, f);
   }
};

// a partial specialization for LHS = RowMajor * ColMajor
// would be rude here - it would force specializations on the
// LHS or the value_type to also follow suit.
// Instead, we forward to AssignProduct2_StrideStrideStride
// and specialize further there.

template <typename LHS, typename LHSOrient, 
	  typename M1, typename M1Orient,
	  typename M2, typename M2Orient,
          typename Nested>
struct AssignProduct2_StrideStrideStride;

template <typename LHS, typename LHSOrient,
	  typename M1, typename M2,
          typename Nested>
struct AssignProduct2_StrideStrideStride<LHS, LHSOrient, M1, RowMajor, M2, ColMajor,
                                         Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested f)
   {
      typename const_iterator<M1>::type I = iterate(m1);
      while (I)
      {
	 typename const_iterator<M2>::type J = iterate(m2);
	 while (J)
	 {
	    set_element(lhs, I.index(), J.index(), parallel_prod(*I, *J, f));
	    ++J;
	 }
	 ++I;
      }
   }
};

template <typename LHS, typename LHSOrient,
	  typename M1, typename M2, typename Nested>
struct AssignProduct2_StrideStrideStride<LHS, LHSOrient, M1, ColMajor, M2, ColMajor,
                                         Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested f)
   {
      assign_product2(lhs, swap_sort_order(m1), m2, f);
   }
};

template <typename LHS, typename LHSOrient,
	  typename M1, typename M2, typename Nested>
struct AssignProduct2_StrideStrideStride<LHS, LHSOrient, M1, RowMajor, M2, RowMajor,
                                         Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested f)
   {
      assign_product2(lhs, m1, swap_sort_order(m2), f);
   }
};

template <typename LHS, typename LHSOrient,
	  typename M1, typename M2, typename Nested>
struct AssignProduct2_StrideStrideStride<LHS, LHSOrient, M1, ColMajor, M2, RowMajor,
                                         Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested f)
   {
      assign_product2(lhs, swap_sort_order(m1), swap_sort_order(m2), f);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1Orient, typename M1i,
	  typename M2v, typename M2Orient, typename M2i>
struct AssignProduct2<LHS, M1, M2, Nested,
		      Concepts::StrideMatrix<LHSv, LHSOrient, LHSi>,
		      Concepts::StrideMatrix<M1v, M1Orient, M1i>,
		      Concepts::StrideMatrix<M2v, M2Orient, M2i>>
: AssignProduct2_StrideStrideStride<LHS, LHSOrient, 
                                    M1, M1Orient, 
                                    M2, M2Orient,
                                    Nested> {};

// assign sparse

template <typename LHS,
	  typename M1, typename M1Orient,
	  typename M2, typename M2Orient,
          typename Nested>
struct AssignProduct2_Sparse;

template <typename LHS,
	  typename M1, typename M2, typename Nested>
struct AssignProduct2_Sparse<LHS, M1, RowMajor, M2, ColMajor, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      zero_all(lhs);
      typename const_iterator<M1>::type I = iterate(m1);
      while (I)
      {
	 typename const_iterator<M2>::type J = iterate(m2);
	 while (J)
	 {
	    add_element_check_if_zero(lhs, I.index(), J.index(), 
                                      parallel_prod(*I, *J, f));
	    ++J;
	 }
	 ++I;
      }
   }
};

template <typename LHS,
	  typename M1, typename M2, typename Nested>
struct AssignProduct2_Sparse<LHS, M1, ColMajor, M2, ColMajor, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      assign_product2(lhs, Temp1, m2, f);
   }
};

template <typename LHS,
	  typename M1, typename M2, typename Nested>
struct AssignProduct2_Sparse<LHS, M1, ColMajor, M2, RowMajor, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      assign_product2(lhs, Temp1, Temp2, f);
   }
};

template <typename LHS,
	  typename M1, typename M2, typename Nested>
struct AssignProduct2_Sparse<LHS, M1, RowMajor, M2, RowMajor, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      assign_product2(lhs, m1, Temp2, f);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1Orient, typename M1i,
	  typename M2v, typename M2Orient, typename M2i>
struct AssignProduct2<LHS, M1, M2, Nested,
		      Concepts::LocalMatrix<LHSv, LHSi>,
		      Concepts::CompressedOuterMatrix<M1v, M1Orient, M1i>,
		      Concepts::CompressedOuterMatrix<M2v, M2Orient, M2i>>
: AssignProduct2_Sparse<LHS, M1, M1Orient, M2, M2Orient, Nested> {};

// assign mixed sparse/dense

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AssignProduct2<LHS, M1, M2, Nested,
		      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
		      Concepts::CompressedOuterMatrix<M1v, RowMajor, M1i>,
		      Concepts::DenseMatrix<M2v, ColMajor, M2i>>
: AssignProduct2_Sparse<LHS, M1, RowMajor, M2, ColMajor, Nested> {};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AssignProduct2<LHS, M1, M2, Nested,
		      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
		      Concepts::CompressedOuterMatrix<M1v, ColMajor, M1i>,
		      Concepts::DenseMatrix<M2v, ColMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      assign_product2(lhs, Temp1, m2, f);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AssignProduct2<LHS, M1, M2, Nested,
		      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
		      Concepts::CompressedOuterMatrix<M1v, ColMajor, M1i>,
		      Concepts::DenseMatrix<M2v, RowMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      assign_product2(lhs, Temp1, swap_sort_order(m2), f);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AssignProduct2<LHS, M1, M2, Nested,
		      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
		      Concepts::CompressedOuterMatrix<M1v, RowMajor, M1i>,
		      Concepts::DenseMatrix<M2v, RowMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      assign_product2(lhs, m1, swap_sort_order(m2), f);
   }
};


template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AssignProduct2<LHS, M1, M2, Nested,
		      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
		      Concepts::DenseMatrix<M1v, RowMajor, M1i>,
		      Concepts::CompressedOuterMatrix<M2v, ColMajor, M2i>>
: AssignProduct2_Sparse<LHS, M1, RowMajor, M2, ColMajor, Nested> {};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AssignProduct2<LHS, M1, M2, Nested,
		      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
		      Concepts::DenseMatrix<M1v, ColMajor, M1i>,
		      Concepts::CompressedOuterMatrix<M2v, ColMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      assign_product2(lhs, swap_sort_order(m1), m2, f);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
					      typename M2v, typename M2i>
					      struct AssignProduct2<LHS, M1, M2, Nested,
					      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
					      Concepts::DenseMatrix<M1v, ColMajor, M1i>,
					      Concepts::CompressedOuterMatrix<M2v, RowMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      assign_product2(lhs, swap_sort_order(m1), Temp2, f);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AssignProduct2<LHS, M1, M2, Nested,
		      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
		      Concepts::DenseMatrix<M1v, RowMajor, M1i>,
		      Concepts::CompressedOuterMatrix<M2v, RowMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      assign_product2(lhs, m1, Temp2, f);
   }
};

//
// add_product2
//

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSInterface = typename interface<LHS>::type,
	  typename M1Interface = typename interface<M1>::type,
	  typename M2Interface = typename interface<M2>::type>
struct AddProduct2 {};

template <typename LHS, typename M1, typename M2, typename Nested>
inline
void add_product2(LHS& lhs, M1 const& m1, M2 const& m2, Nested f)
{
   AddProduct2<LHS, M1, M2, Nested>::apply(lhs, m1, m2, f);
}

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AddProduct2<LHS, M1, M2, Nested,
		      Concepts::MatrixExpression<LHSv, LHSi>,
		      Concepts::MatrixExpression<M1v, M1i>,
		      Concepts::MatrixExpression<M2v, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested f)
   {
      add_product2(lhs, eval_expression(m1), eval_expression(m2), f);
   }
};

// add dense

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1Orient, typename M1i,
	  typename M2v, typename M2Orient, typename M2i>
struct AddProduct2<LHS, M1, M2, Nested,
		      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
		      Concepts::DenseMatrix<M1v, M1Orient, M1i>,
		      Concepts::DenseMatrix<M2v, M2Orient, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested f)
   {
      Matrix<M1v, RowMajor> Temp1(m1);
      Matrix<M2v, ColMajor> Temp2(m2);
      Matrix<LHSv, RowMajor> TempL(size1(m1), size2(m2));
      assign_product2(TempL, Temp1, Temp2, f);
      add(lhs, TempL);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1Orient, typename M1i,
	  typename M2v, typename M2Orient, typename M2i>
struct AddProduct2<LHS, M1, M2, Nested,
		      Concepts::StrideMatrix<LHSv, LHSOrient, LHSi>,
		      Concepts::DenseMatrix<M1v, M1Orient, M1i>,
		      Concepts::DenseMatrix<M2v, M2Orient, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested f)
   {
      Matrix<M1v, RowMajor> Temp1(m1);
      Matrix<M2v, ColMajor> Temp2(m2);
      add_product2(lhs, Temp1, Temp2, f);
   }
};

// a partial specialization for LHS = RowMajor * ColMajor
// would be rude here - it would force specializations on the
// LHS or the value_type to also follow suit.
// Instead, we forward to AddProduct2_StrideStrideStride
// and specialize further there.

template <typename LHS, typename LHSOrient, 
	  typename M1, typename M1Orient,
	  typename M2, typename M2Orient,
          typename Nested>
struct AddProduct2_StrideStrideStride;

template <typename LHS, typename LHSOrient,
	  typename M1, typename M2,
          typename Nested>
struct AddProduct2_StrideStrideStride<LHS, LHSOrient, M1, RowMajor, M2, ColMajor,
                                      Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested f)
   {
      typename const_iterator<M1>::type I = iterate(m1);
      while (I)
      {
	 typename const_iterator<M2>::type J = iterate(m2);
	 while (J)
	 {
	    add_element(lhs, I.index(), J.index(), parallel_prod(*I, *J, f));
	    ++J;
	 }
	 ++I;
      }
   }
};

template <typename LHS, typename LHSOrient,
	  typename M1, typename M2, typename Nested>
struct AddProduct2_StrideStrideStride<LHS, LHSOrient, M1, ColMajor, M2, ColMajor,
                                      Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested f)
   {
      add_product2(lhs, swap_sort_order(m1), m2, f);
   }
};

template <typename LHS, typename LHSOrient,
	  typename M1, typename M2, typename Nested>
struct AddProduct2_StrideStrideStride<LHS, LHSOrient, M1, RowMajor, M2, RowMajor,
                                      Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested f)
   {
      add_product2(lhs, m1, swap_sort_order(m2), f);
   }
};

template <typename LHS, typename LHSOrient,
	  typename M1, typename M2, typename Nested>
struct AddProduct2_StrideStrideStride<LHS, LHSOrient, M1, ColMajor, M2, RowMajor,
                                      Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested f)
   {
      add_product2(lhs, swap_sort_order(m1), swap_sort_order(m2), f);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1Orient, typename M1i,
	  typename M2v, typename M2Orient, typename M2i>
struct AddProduct2<LHS, M1, M2, Nested,
		      Concepts::StrideMatrix<LHSv, LHSOrient, LHSi>,
		      Concepts::StrideMatrix<M1v, M1Orient, M1i>,
		      Concepts::StrideMatrix<M2v, M2Orient, M2i>>
: AddProduct2_StrideStrideStride<LHS, LHSOrient, M1, M1Orient, M2, M2Orient,
                                 Nested> {};

// add sparse

template <typename LHS,
	  typename M1, typename M1Orient,
	  typename M2, typename M2Orient, typename Nested>
struct AddProduct2_Sparse;

template <typename LHS,
	  typename M1, typename M2, typename Nested>
struct AddProduct2_Sparse<LHS, M1, RowMajor, M2, ColMajor, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      typename const_iterator<M1>::type I = iterate(m1);
      while (I)
      {
	 typename const_iterator<M2>::type J = iterate(m2);
	 while (J)
	 {
	    add_element_check_if_zero(lhs, I.index(), J.index(), 
                                      parallel_prod(*I, *J, f));
	    ++J;
	 }
	 ++I;
      }
   }
};

template <typename LHS,
	  typename M1, typename M2, typename Nested>
struct AddProduct2_Sparse<LHS, M1, ColMajor, M2, ColMajor, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      add_product2(lhs, Temp1, m2, f);
   }
};

template <typename LHS,
	  typename M1, typename M2, typename Nested>
struct AddProduct2_Sparse<LHS, M1, ColMajor, M2, RowMajor, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      add_product2(lhs, Temp1, Temp2, f);
   }
};

template <typename LHS,
	  typename M1, typename M2, typename Nested>
struct AddProduct2_Sparse<LHS, M1, RowMajor, M2, RowMajor, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      add_product2(lhs, m1, Temp2, f);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1Orient, typename M1i,
	  typename M2v, typename M2Orient, typename M2i>
struct AddProduct2<LHS, M1, M2, Nested,
		   Concepts::LocalMatrix<LHSv, LHSi>,
		   Concepts::CompressedOuterMatrix<M1v, M1Orient, M1i>,
		   Concepts::CompressedOuterMatrix<M2v, M2Orient, M2i>>
: AddProduct2_Sparse<LHS, M1, M1Orient, M2, M2Orient, Nested> {};

// add mixed sparse/dense

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AddProduct2<LHS, M1, M2, Nested,
		      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
		      Concepts::CompressedOuterMatrix<M1v, RowMajor, M1i>,
		      Concepts::DenseMatrix<M2v, ColMajor, M2i>>
: AddProduct2_Sparse<LHS, M1, RowMajor, M2, ColMajor, Nested> {};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AddProduct2<LHS, M1, M2, Nested,
		      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
		      Concepts::CompressedOuterMatrix<M1v, ColMajor, M1i>,
		      Concepts::DenseMatrix<M2v, ColMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      add_product2(lhs, Temp1, m2, f);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AddProduct2<LHS, M1, M2, Nested,
		      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
		      Concepts::CompressedOuterMatrix<M1v, ColMajor, M1i>,
		      Concepts::DenseMatrix<M2v, RowMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      add_product2(lhs, Temp1, swap_sort_order(m2), f);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AddProduct2<LHS, M1, M2, Nested,
		      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
		      Concepts::CompressedOuterMatrix<M1v, RowMajor, M1i>,
		      Concepts::DenseMatrix<M2v, RowMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      add_product2(lhs, m1, swap_sort_order(m2), f);
   }
};


template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AddProduct2<LHS, M1, M2, Nested,
		      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
		      Concepts::DenseMatrix<M1v, RowMajor, M1i>,
		      Concepts::CompressedOuterMatrix<M2v, ColMajor, M2i>>
: AddProduct2_Sparse<LHS, M1, RowMajor, M2, ColMajor, Nested> {};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AddProduct2<LHS, M1, M2, Nested,
		      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
		      Concepts::DenseMatrix<M1v, ColMajor, M1i>,
		      Concepts::CompressedOuterMatrix<M2v, ColMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      add_product2(lhs, swap_sort_order(m1), m2, f);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AddProduct2<LHS, M1, M2, Nested,
		      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
		      Concepts::DenseMatrix<M1v, ColMajor, M1i>,
		      Concepts::CompressedOuterMatrix<M2v, RowMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      add_product2(lhs, swap_sort_order(m1), Temp2, f);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AddProduct2<LHS, M1, M2, Nested,
		      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
		      Concepts::DenseMatrix<M1v, RowMajor, M1i>,
		      Concepts::CompressedOuterMatrix<M2v, RowMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      add_product2(lhs, m1, Temp2, f);
   }
};

//
// subtract_product2
//

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSInterface = typename interface<LHS>::type,
	  typename M1Interface = typename interface<M1>::type,
	  typename M2Interface = typename interface<M2>::type>
struct SubtractProduct2 {};

template <typename LHS, typename M1, typename M2, typename Nested>
inline
void subtract_product2(LHS& lhs, M1 const& m1, M2 const& m2, Nested f)
{
   SubtractProduct2<LHS, M1, M2, Nested>::apply(lhs, m1, m2, f);
}

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct SubtractProduct2<LHS, M1, M2, Nested,
		      Concepts::MatrixExpression<LHSv, LHSi>,
		      Concepts::MatrixExpression<M1v, M1i>,
		      Concepts::MatrixExpression<M2v, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested f)
   {
      subtract_product2(lhs, eval_expression(m1), eval_expression(m2), f);
   }
};

// subtract dense

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1Orient, typename M1i,
	  typename M2v, typename M2Orient, typename M2i>
struct SubtractProduct2<LHS, M1, M2, Nested,
		      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
		      Concepts::DenseMatrix<M1v, M1Orient, M1i>,
		      Concepts::DenseMatrix<M2v, M2Orient, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested f)
   {
      Matrix<M1v, RowMajor> Temp1(m1);
      Matrix<M2v, ColMajor> Temp2(m2);
      Matrix<LHSv, RowMajor> TempL(size1(m1), size2(m2));
      assign_product2(TempL, Temp1, Temp2, f);
      subtract(lhs, TempL);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1Orient, typename M1i,
	  typename M2v, typename M2Orient, typename M2i>
struct SubtractProduct2<LHS, M1, M2, Nested,
		      Concepts::StrideMatrix<LHSv, LHSOrient, LHSi>,
		      Concepts::DenseMatrix<M1v, M1Orient, M1i>,
		      Concepts::DenseMatrix<M2v, M2Orient, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested f)
   {
      Matrix<M1v, RowMajor> Temp1(m1);
      Matrix<M2v, ColMajor> Temp2(m2);
      subtract_product2(lhs, Temp1, Temp2, f);
   }
};

// a partial specialization for LHS = RowMajor * ColMajor
// would be rude here - it would force specializations on the
// LHS or the value_type to also follow suit.
// Instead, we forward to SubtractProduct2_StrideStrideStride
// and specialize further there.

template <typename LHS, typename LHSOrient, 
	  typename M1, typename M1Orient,
	  typename M2, typename M2Orient,
          typename Nested>
struct SubtractProduct2_StrideStrideStride;

template <typename LHS, typename LHSOrient,
	  typename M1, typename M2,
          typename Nested>
struct SubtractProduct2_StrideStrideStride<LHS, LHSOrient, M1, RowMajor, M2, ColMajor,
                                      Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested f)
   {
      typename const_iterator<M1>::type I = iterate(m1);
      while (I)
      {
	 typename const_iterator<M2>::type J = iterate(m2);
	 while (J)
	 {
	    subtract_element(lhs, I.index(), J.index(), parallel_prod(*I, *J, f));
	    ++J;
	 }
	 ++I;
      }
   }
};

template <typename LHS, typename LHSOrient,
	  typename M1, typename M2, typename Nested>
struct SubtractProduct2_StrideStrideStride<LHS, LHSOrient, M1, ColMajor, M2, ColMajor,
                                      Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested f)
   {
      subtract_product2(lhs, swap_sort_order(m1), m2, f);
   }
};

template <typename LHS, typename LHSOrient,
	  typename M1, typename M2, typename Nested>
struct SubtractProduct2_StrideStrideStride<LHS, LHSOrient, M1, RowMajor, M2, RowMajor,
                                      Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested f)
   {
      subtract_product2(lhs, m1, swap_sort_order(m2), f);
   }
};

template <typename LHS, typename LHSOrient,
	  typename M1, typename M2, typename Nested>
struct SubtractProduct2_StrideStrideStride<LHS, LHSOrient, M1, ColMajor, M2, RowMajor,
                                      Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested f)
   {
      subtract_product2(lhs, swap_sort_order(m1), swap_sort_order(m2), f);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1Orient, typename M1i,
	  typename M2v, typename M2Orient, typename M2i>
struct SubtractProduct2<LHS, M1, M2, Nested,
			Concepts::StrideMatrix<LHSv, LHSOrient, LHSi>,
			Concepts::StrideMatrix<M1v, M1Orient, M1i>,
			Concepts::StrideMatrix<M2v, M2Orient, M2i>>
: SubtractProduct2_StrideStrideStride<LHS, LHSOrient, M1, M1Orient, M2, M2Orient,
                                 Nested> {};

// subtract sparse

template <typename LHS,
	  typename M1, typename M1Orient,
	  typename M2, typename M2Orient, typename Nested>
struct SubtractProduct2_Sparse;

template <typename LHS,
	  typename M1, typename M2, typename Nested>
struct SubtractProduct2_Sparse<LHS, M1, RowMajor, M2, ColMajor, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      typename const_iterator<M1>::type I = iterate(m1);
      while (I)
      {
	 typename const_iterator<M2>::type J = iterate(m2);
	 while (J)
	 {
	    subtract_element_check_if_zero(lhs, I.index(), J.index(), 
                                      parallel_prod(*I, *J, f));
	    ++J;
	 }
	 ++I;
      }
   }
};

template <typename LHS,
	  typename M1, typename M2, typename Nested>
struct SubtractProduct2_Sparse<LHS, M1, ColMajor, M2, ColMajor, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      subtract_product2(lhs, Temp1, m2, f);
   }
};

template <typename LHS,
	  typename M1, typename M2, typename Nested>
struct SubtractProduct2_Sparse<LHS, M1, ColMajor, M2, RowMajor, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      subtract_product2(lhs, Temp1, Temp2, f);
   }
};

template <typename LHS,
	  typename M1, typename M2, typename Nested>
struct SubtractProduct2_Sparse<LHS, M1, RowMajor, M2, RowMajor, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      subtract_product2(lhs, m1, Temp2, f);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1Orient, typename M1i,
	  typename M2v, typename M2Orient, typename M2i>
struct SubtractProduct2<LHS, M1, M2, Nested,
			Concepts::LocalMatrix<LHSv, LHSi>,
			Concepts::CompressedOuterMatrix<M1v, M1Orient, M1i>,
			Concepts::CompressedOuterMatrix<M2v, M2Orient, M2i>>
: SubtractProduct2_Sparse<LHS, M1, M1Orient, M2, M2Orient, Nested> {};

// subtract mixed sparse/dense

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct SubtractProduct2<LHS, M1, M2, Nested,
		      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
		      Concepts::CompressedOuterMatrix<M1v, RowMajor, M1i>,
		      Concepts::DenseMatrix<M2v, ColMajor, M2i>>
: SubtractProduct2_Sparse<LHS, M1, RowMajor, M2, ColMajor, Nested> {};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct SubtractProduct2<LHS, M1, M2, Nested,
		      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
		      Concepts::CompressedOuterMatrix<M1v, ColMajor, M1i>,
		      Concepts::DenseMatrix<M2v, ColMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      subtract_product2(lhs, Temp1, m2, f);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct SubtractProduct2<LHS, M1, M2, Nested,
		      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
		      Concepts::CompressedOuterMatrix<M1v, ColMajor, M1i>,
		      Concepts::DenseMatrix<M2v, RowMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      subtract_product2(lhs, Temp1, swap_sort_order(m2), f);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct SubtractProduct2<LHS, M1, M2, Nested,
		      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
		      Concepts::CompressedOuterMatrix<M1v, RowMajor, M1i>,
		      Concepts::DenseMatrix<M2v, RowMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      subtract_product2(lhs, m1, swap_sort_order(m2), f);
   }
};


template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct SubtractProduct2<LHS, M1, M2, Nested,
		      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
		      Concepts::DenseMatrix<M1v, RowMajor, M1i>,
		      Concepts::CompressedOuterMatrix<M2v, ColMajor, M2i>>
: SubtractProduct2_Sparse<LHS, M1, RowMajor, M2, ColMajor, Nested> {};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct SubtractProduct2<LHS, M1, M2, Nested,
		      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
		      Concepts::DenseMatrix<M1v, ColMajor, M1i>,
		      Concepts::CompressedOuterMatrix<M2v, ColMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      subtract_product2(lhs, swap_sort_order(m1), m2, f);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct SubtractProduct2<LHS, M1, M2, Nested,
		      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
		      Concepts::DenseMatrix<M1v, ColMajor, M1i>,
		      Concepts::CompressedOuterMatrix<M2v, RowMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      subtract_product2(lhs, swap_sort_order(m1), Temp2, f);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct SubtractProduct2<LHS, M1, M2, Nested,
		      Concepts::DenseMatrix<LHSv, LHSOrient, LHSi>,
		      Concepts::DenseMatrix<M1v, RowMajor, M1i>,
		      Concepts::CompressedOuterMatrix<M2v, RowMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      subtract_product2(lhs, m1, Temp2, f);
   }
};

} // namespace LinearAlgebra

#if !defined(LINEARALGEBRA_NO_BLAS)
#include "matrixproductblas.h"
#endif

#endif
