// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/matrixproductdiagonal.h
//
// Copyright (C) 2009-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
   Created 2009-03-17 by Ian McCulloch

   Specializations of the matrix-matrix primitives involving diagonal matrices

*/

#if !defined(MPTOOLKIT_LINEARALGEBRA_MATRIXPRODUCTDIAGONAL_H)
#define MPTOOLKIT_LINEARALGEBRA_MATRIXPRODUCTDIAGONAL_H

#include "matrixproductoperations.h"
#include "diagonalmatrix.h"

namespace LinearAlgebra
{

//
// assign_product2
//

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AssignProduct2<LHS, M1, M2, Nested,
		      Concepts::StrideMatrix<LHSv, ColMajor, LHSi>,
		      Concepts::DiagonalMatrix<M1v, M1i>,
		      Concepts::StrideMatrix<M2v, ColMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      DEBUG_PRECONDITION_EQUAL(size2(m1), size1(m2));
      DEBUG_PRECONDITION_EQUAL(size1(lhs), size1(m1));
      DEBUG_PRECONDITION_EQUAL(size2(lhs), size2(m2));
      typename const_iterator<M2>::type J = iterate(m2);
      typename iterator<LHS>::type K = iterate(lhs);
      while (J)
      {
	 *K = transform(m1.diagonal(), *J, f);
	 ++J; ++K;
      }
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AssignProduct2<LHS, M1, M2, Nested,
		      Concepts::StrideMatrix<LHSv, RowMajor, LHSi>,
		      Concepts::DiagonalMatrix<M1v, M1i>,
		      Concepts::StrideMatrix<M2v, ColMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      assign_product2(swap_sort_order(lhs), m1, m2, f);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AssignProduct2<LHS, M1, M2, Nested,
		      Concepts::StrideMatrix<LHSv, RowMajor, LHSi>,
		      Concepts::DiagonalMatrix<M1v, M1i>,
		      Concepts::StrideMatrix<M2v, RowMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      assign_product2(swap_sort_order(lhs), m1, swap_sort_order(m2), f);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AssignProduct2<LHS, M1, M2, Nested,
		      Concepts::StrideMatrix<LHSv, ColMajor, LHSi>,
		      Concepts::DiagonalMatrix<M1v, M1i>,
		      Concepts::StrideMatrix<M2v, RowMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      assign_product2(lhs, m1, swap_sort_order(m2), f);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AssignProduct2<LHS, M1, M2, Nested,
		      Concepts::StrideMatrix<LHSv, RowMajor, LHSi>,
		      Concepts::StrideMatrix<M1v, RowMajor, M1i>,
		      Concepts::DiagonalMatrix<M2v, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      DEBUG_PRECONDITION_EQUAL(size2(m1), size1(m2));
      DEBUG_PRECONDITION_EQUAL(size1(lhs), size1(m1));
      DEBUG_PRECONDITION_EQUAL(size2(lhs), size2(m2));
      typename const_iterator<M1>::type I = iterate(m1);
      typename iterator<LHS>::type K = iterate(lhs);
      while (I)
      {
	 *K = transform(*I, m2.diagonal(), f);
	 ++I; ++K;
      }
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AssignProduct2<LHS, M1, M2, Nested,
		      Concepts::StrideMatrix<LHSv, RowMajor, LHSi>,
		      Concepts::StrideMatrix<M1v, ColMajor, M1i>,
		      Concepts::DiagonalMatrix<M2v, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      assign_product2(lhs, swap_sort_order(m1), m2, f);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AssignProduct2<LHS, M1, M2, Nested,
		      Concepts::StrideMatrix<LHSv, ColMajor, LHSi>,
		      Concepts::StrideMatrix<M1v, ColMajor, M1i>,
		      Concepts::DiagonalMatrix<M2v, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      assign_product2(swap_sort_order(lhs), swap_sort_order(m1), m2, f);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AssignProduct2<LHS, M1, M2, Nested,
		      Concepts::StrideMatrix<LHSv, ColMajor, LHSi>,
		      Concepts::StrideMatrix<M1v, RowMajor, M1i>,
		      Concepts::DiagonalMatrix<M2v, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      assign_product2(swap_sort_order(lhs), m1, m2, f);
   }
};

//
// add_product2
//

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AddProduct2<LHS, M1, M2, Nested,
		      Concepts::StrideMatrix<LHSv, ColMajor, LHSi>,
		      Concepts::DiagonalMatrix<M1v, M1i>,
		      Concepts::StrideMatrix<M2v, ColMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      DEBUG_PRECONDITION_EQUAL(size2(m1), size1(m2));
      DEBUG_PRECONDITION_EQUAL(size1(lhs), size1(m1));
      DEBUG_PRECONDITION_EQUAL(size2(lhs), size2(m2));
      typename IterateDiagonal<M1>::result_type I = iterate(m1.diagonal());
      typename const_iterator<M2>::type J = iterate(m2);
      typename iterator<LHS>::type K = iterate(lhs);
      while (I)
      {
	 *K += transform(*I, *J, f);
	 ++I; ++J; ++K;
      }
      //SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      //add_product2(lhs, Temp1, swap_sort_order(m2), f);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AddProduct2<LHS, M1, M2, Nested,
		      Concepts::StrideMatrix<LHSv, RowMajor, LHSi>,
		      Concepts::DiagonalMatrix<M1v, M1i>,
		      Concepts::StrideMatrix<M2v, ColMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      add_product2(swap_sort_order(lhs), m1, m2, f);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AddProduct2<LHS, M1, M2, Nested,
		      Concepts::StrideMatrix<LHSv, RowMajor, LHSi>,
		      Concepts::DiagonalMatrix<M1v, M1i>,
		      Concepts::StrideMatrix<M2v, RowMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      add_product2(swap_sort_order(lhs), m1, swap_sort_order(m2), f);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AddProduct2<LHS, M1, M2, Nested,
		      Concepts::StrideMatrix<LHSv, ColMajor, LHSi>,
		      Concepts::DiagonalMatrix<M1v, M1i>,
		      Concepts::StrideMatrix<M2v, RowMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      add_product2(lhs, m1, swap_sort_order(m2), f);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AddProduct2<LHS, M1, M2, Nested,
		   Concepts::StrideMatrix<LHSv, RowMajor, LHSi>,
		   Concepts::StrideMatrix<M1v, RowMajor, M1i>,
		   Concepts::DiagonalMatrix<M2v, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      DEBUG_PRECONDITION_EQUAL(size2(m1), size1(m2));
      DEBUG_PRECONDITION_EQUAL(size1(lhs), size1(m1));
      DEBUG_PRECONDITION_EQUAL(size2(lhs), size2(m2));
      typename const_iterator<M1>::type I = iterate(m2);
      typename iterator<LHS>::type K = iterate(lhs);
      while (I)
      {
	 *K += transform(*I, m2.diagonal(), f);
	 ++I; ++K;
      }
      //SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      //add_product2(lhs, Temp1, swap_sort_order(m2), f);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AddProduct2<LHS, M1, M2, Nested,
		   Concepts::StrideMatrix<LHSv, RowMajor, LHSi>,
		   Concepts::StrideMatrix<M1v, ColMajor, M1i>,
		   Concepts::DiagonalMatrix<M2v, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      add_product2(lhs, m1, swap_sort_order(m2), f);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AddProduct2<LHS, M1, M2, Nested,
		   Concepts::StrideMatrix<LHSv, ColMajor, LHSi>,
		   Concepts::StrideMatrix<M1v, ColMajor, M1i>,
		   Concepts::DiagonalMatrix<M2v, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      add_product2(swap_sort_order(lhs), m1, swap_sort_order(m2), f);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AddProduct2<LHS, M1, M2, Nested,
		   Concepts::StrideMatrix<LHSv, ColMajor, LHSi>,
		   Concepts::StrideMatrix<M1v, RowMajor, M1i>,
		   Concepts::DiagonalMatrix<M2v, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      add_product2(swap_sort_order(lhs), m1, m2, f);
   }
};

//
// subtract_product2
//

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct SubtractProduct2<LHS, M1, M2, Nested,
			Concepts::StrideMatrix<LHSv, ColMajor, LHSi>,
			Concepts::DiagonalMatrix<M1v, M1i>,
			Concepts::StrideMatrix<M2v, ColMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      DEBUG_PRECONDITION_EQUAL(size2(m1), size1(m2));
      DEBUG_PRECONDITION_EQUAL(size1(lhs), size1(m1));
      DEBUG_PRECONDITION_EQUAL(size2(lhs), size2(m2));
      typename const_iterator<M2>::type J = iterate(m2);
      typename iterator<LHS>::type K = iterate(lhs);
      while (J)
      {
	 *K -= transform(m1.diagonal(), *J, f);
	 ++J; ++K;
      }
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct SubtractProduct2<LHS, M1, M2, Nested,
			Concepts::StrideMatrix<LHSv, RowMajor, LHSi>,
			Concepts::DiagonalMatrix<M1v, M1i>,
			Concepts::StrideMatrix<M2v, ColMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      subtract_product2(swap_sort_order(lhs), m1, m2, f);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct SubtractProduct2<LHS, M1, M2, Nested,
			Concepts::StrideMatrix<LHSv, RowMajor, LHSi>,
			Concepts::DiagonalMatrix<M1v, M1i>,
			Concepts::StrideMatrix<M2v, RowMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      subtract_product2(swap_sort_order(lhs), m1, swap_sort_order(m2), f);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct SubtractProduct2<LHS, M1, M2, Nested,
			Concepts::StrideMatrix<LHSv, ColMajor, LHSi>,
			Concepts::DiagonalMatrix<M1v, M1i>,
			Concepts::StrideMatrix<M2v, RowMajor, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      subtract_product2(lhs, m1, swap_sort_order(m2), f);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct SubtractProduct2<LHS, M1, M2, Nested,
			Concepts::StrideMatrix<LHSv, RowMajor, LHSi>,
			Concepts::StrideMatrix<M1v, RowMajor, M1i>,
			Concepts::DiagonalMatrix<M2v, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      DEBUG_PRECONDITION_EQUAL(size2(m1), size1(m2));
      DEBUG_PRECONDITION_EQUAL(size1(lhs), size1(m1));
      DEBUG_PRECONDITION_EQUAL(size2(lhs), size2(m2));
      typename const_iterator<M1>::type I = iterate(m1);
      typename iterator<LHS>::type K = iterate(lhs);
      while (I)
      {
	 *K -= transform(*I, m2.diagonal(), f);
	 ++I; ++K;
      }
      //SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      //subtract_product2(lhs, Temp1, swap_sort_order(m2), f);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct SubtractProduct2<LHS, M1, M2, Nested,
			Concepts::StrideMatrix<LHSv, RowMajor, LHSi>,
			Concepts::StrideMatrix<M1v, ColMajor, M1i>,
			Concepts::DiagonalMatrix<M2v, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      subtract_product2(lhs, m1, swap_sort_order(m2), f);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct SubtractProduct2<LHS, M1, M2, Nested,
			Concepts::StrideMatrix<LHSv, ColMajor, LHSi>,
			Concepts::StrideMatrix<M1v, ColMajor, M1i>,
			Concepts::DiagonalMatrix<M2v, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      subtract_product2(swap_sort_order(lhs), m1, swap_sort_order(m2), f);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct SubtractProduct2<LHS, M1, M2, Nested,
			Concepts::StrideMatrix<LHSv, ColMajor, LHSi>,
			Concepts::StrideMatrix<M1v, RowMajor, M1i>,
			Concepts::DiagonalMatrix<M2v, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      subtract_product2(swap_sort_order(lhs), m1, m2, f);
   }
};


} // namespace LinearAlgebra

#endif
