// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/matrixinterface.h
//
// Copyright (C) 2005-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
  matrixinterface.h

  New version of the interface classes for the matrix library.

  Created 2005-01-14 Ian McCulloch
*/

#if !defined(MATRIXINTERFACE_H_KJHIU4YR873Y68YFO87YRO8743)
#define MATRIXINTERFACE_H_KJHIU4YR873Y68YFO87YRO8743

#include "vectorinterface.h"

namespace LinearAlgebra
{

// ContiguousRows means that the row vectors are contiguous (ie. stride2() == 1),
// but stride1() != size2() [ie. there may be a gap between the end of one
// row and the start of the next].
//
// ContiguousRowMajor means that the storage is completely contiguous,
// ie. stride2() == 1 && stride1() == size2().
//
// For a coordinate matrix, the nested interface must be Compressed (and refinements).
// Maybe we would want to add variants on Ordered?  Eg, RowMajorOrdering,
// ColMajorOrdering ?  overlap with XXXMajorMatrix and rvalue vs lvalue issues here.
//
// Possibly, distributed matrices would fit into RowMajorMatrix/ColMajorMatrix
// with a distributed nested interface; not sure how this would work yet.
//
// Banded and Dense matrices probably should have fixed size variants.
//
// A triangular matrix, if we added one, would inherit from CompressedMatrix.
//
// A symmetric matrix could occur in sparse or dense variations; possibly
// have a DenseSymmetric inherit from DenseMatrix, and a SparseSymmetric
// inherit from CompressedMatrix?
//
// TODO: specify iterator sort order.
//
// TODO: specify constraints on nested interfaces.

//
// 2015-10-03: The InjectiveMatrix isn't used - decided to make
// SparseMatrix injective.

// tags for row/column major
struct RowMajor { typedef RowMajor type; };
struct ColMajor { typedef ColMajor type; };

//
// SwapOrientation - metafunction for interchanging RowMajor <-> ColMajor
//

template <typename Orientation>
struct SwapOrientation;

// this depends on ColMajor having a nested typedef type = ColMajor
template <> struct SwapOrientation<RowMajor> : ColMajor {};
template <> struct SwapOrientation<ColMajor> : RowMajor {};

namespace Concepts
{

// tags

template <typename T>
struct ExpressionTag {};

template <typename T>
struct LocalTag {};

template <typename T>
struct InjectiveTag {};

template <typename T>
struct SparseTag {};

template <typename T>
struct CoordinateTag {};

template <typename T>
struct DiagonalTag {};

template <typename T>
struct ScalarTag {};

template <typename Orientation, typename T>
struct CompressedOuterTag {};

template <typename Orientation, typename T>
struct DenseTag {};

template <typename Orientation, typename T>
struct StrideTag {};

template <typename Orientation, typename T>
struct ContiguousTag {};

// concept heirachy

// top level concept
template <typename Value, typename T>
struct AnyMatrix {};

template <typename Value, typename T>
using MatrixExpression = AnyMatrix<Value, ExpressionTag<T>>;

template <typename Value, typename T>
using LocalMatrix = MatrixExpression<Value, LocalTag<T>>;

template <typename Value, typename T>
using InjectiveMatrix = LocalMatrix<Value, InjectiveTag<T>>;

template <typename Value, typename T>
using SparseMatrix = InjectiveMatrix<Value, SparseTag<T>>;

template <typename T, typename Value>
using DiagonalMatrix = SparseMatrix<T, DiagonalTag<Value>>;

template <typename T, typename Value>
using ScalarMatrix = DiagonalMatrix<T, ScalarTag<Value>>;

template <typename Value, typename T>
using CoordinateMatrix = SparseMatrix<Value, CoordinateTag<T>>;

template <typename Value, typename Orientation, typename T>
using CompressedOuterMatrix = SparseMatrix<Value, CompressedOuterTag<Orientation, T>>;

template <typename Value, typename T>
using CompressedRowMajorMatrix = SparseMatrix<Value, CompressedOuterTag<RowMajor, T>>;

template <typename Value, typename T>
using CompressedColMajorMatrix = SparseMatrix<Value, CompressedOuterTag<ColMajor, T>>;

// DenseMatrix

template <typename Value, typename Orientation, typename T>
using DenseMatrix = SparseMatrix<Value, DenseTag<Orientation, T>>;

template <typename Value, typename T>
using DenseRowMajorMatrix = DenseMatrix<Value, RowMajor, T>;

template <typename Value, typename T>
using DenseColMajorMatrix = DenseMatrix<Value, ColMajor, T>;

// StrideMatrix

template <typename Value, typename Orientation, typename T>
using StrideMatrix = DenseMatrix<Value, Orientation, StrideTag<Orientation, T>>;

template <typename Value, typename T>
using StrideRowMajorMatrix = StrideMatrix<Value, RowMajor, T>;

template <typename Value, typename T>
using StrideColMajorMatrix = StrideMatrix<Value, ColMajor, T>;

// ContiguousMatrix

template <typename Value, typename Orientation, typename T>
using ContiguousMatrix = StrideMatrix<Value, Orientation, ContiguousTag<Orientation, T>>;

template <typename Value, typename T>
using ContiguousRowMajorMatrix = ContiguousMatrix<Value, RowMajor, T>;

template <typename Value, typename T>
using ContiguousColMajorMatrix = ContiguousMatrix<Value, ColMajor, T>;

} // namespace Concepts

// is_matrix, boolean function to determine if a type
// has a matrix interface.

namespace Private
{
template <typename T>
struct is_matrix_helper : boost::mpl::false_ {};

template <typename S, typename T>
struct is_matrix_helper<Concepts::AnyMatrix<S,T>> : boost::mpl::true_ {};

template <typename T>
struct is_dense_matrix_helper : boost::mpl::false_ {};

template <typename S, typename T, typename U>
struct is_dense_matrix_helper<Concepts::DenseMatrix<S,T,U>> : boost::mpl::true_ {};

} // namespace Private

template <typename T, typename Enable = void>
struct is_matrix : boost::mpl::false_ {};

template <typename T>
struct is_matrix<T, typename boost::enable_if<exists<typename interface<T>::type> >::type>
   : Private::is_matrix_helper<typename interface<T>::type> {};

template <typename T, typename Enable = void>
struct is_dense_matrix : boost::mpl::false_ {};

template <typename T>
struct is_dense_matrix<T,
   typename boost::enable_if<exists<typename interface<T>::type> >::type>
   : Private::is_dense_matrix_helper<typename interface<T>::type> {};


//
// when allocating temporaries, it is useful to know
// simply whether a type should be sparse or dense.
// this needs a lot of work - should encapsulate both shape & sparseness.
//

struct matrix_abstract_dense
{
   typedef matrix_abstract_dense type;
};

struct matrix_abstract_sparse
{
   typedef matrix_abstract_sparse type;
};

// union of matrices; dense+anything -> dense;  sparse+sparse -> sparse.
template <typename T, typename U>
struct matrix_abstract_or
   : boost::mpl::if_<boost::is_same<typename T::type, typename U::type>, typename T::type,
                     matrix_abstract_dense> {};

// intersection matrices; sparse+anything -> sparse;  dense+dense -> dense.
template <typename T, typename U>
struct matrix_abstract_and
   : boost::mpl::if_<boost::is_same<typename T::type, typename U::type>,
                     typename T::type, matrix_abstract_sparse> {};

// The abstract interfaces for concrete types
template <typename T, typename Tv, typename Orient, typename Ti>
struct abstract_interface_interface<T, Concepts::DenseMatrix<Tv, Orient, Ti>>
   : matrix_abstract_dense {};

template <typename T, typename Tv, typename Ti>
struct abstract_interface_interface<T, Concepts::SparseMatrix<Tv, Ti>>
   : matrix_abstract_sparse {};

//
// iterator categories
//

// generic sparse unordered iterator
struct matrix_iterator_sparse
{
   typedef matrix_iterator_sparse type;
};

// unordered but distinct indices
struct matrix_iterator_injective : matrix_iterator_sparse
{
   typedef matrix_iterator_injective type;
};

// supports lookup via operator()(row,col)  - not yet implemented
struct matrix_iterator_lookup : matrix_iterator_injective
{
   typedef matrix_iterator_lookup type;
};

// dense indices - this should inherit from matrix_iterator_lookup eventually
struct matrix_iterator_dense : matrix_iterator_injective
{
   typedef matrix_iterator_dense type;
};


//
// get_matrix_category
// metafunction for computing the matrix iterator category
// from outer/inner vector iterators.
//

namespace Private
{

typedef char IsSparse[1];
typedef char IsInjective[2];
typedef char IsLookup[3];
typedef char IsDense[4];

template <int Size>
struct category_from_size;

template <> struct category_from_size<1> : matrix_iterator_sparse {};
template <> struct category_from_size<2> : matrix_iterator_injective {};
template <> struct category_from_size<3> : matrix_iterator_lookup {};
template <> struct category_from_size<4> : matrix_iterator_dense {};

template <typename T>
T* MakePtr();   // dummy function, not implemented

IsSparse& TestCategory(...);

IsInjective& TestCategory(vector_iterator_injective const*,
                          vector_iterator_injective const*);

IsDense& TestCategory(vector_iterator_dense const*,
                      vector_iterator_dense const*);

} // namespace Private

template <typename T, typename U>
struct get_matrix_category
   : Private::category_from_size<
        sizeof(Private::TestCategory(Private::MakePtr<T>(), Private::MakePtr<U>()))
     > {};

//
// inner_iterator
// helper metafunction to get the inner iterator type of a matrix
//

template <typename T>
struct inner_iterator : public iterator<typename iterator<T>::type> {};

template <typename T>
struct const_inner_iterator : public iterator<typename const_iterator<T>::type> {};


// some utility functions

template <typename T>
inline
bool is_blas_matrix(T const& M)
{
   return stride1(M) == 1 || stride2(M) == 1;
}

template <typename T>
inline
char blas_trans_col(T const& M)
{
   // note: it is important that we are consistent with comparing stride2() first
   // in blas_trans_col, blas_trans_row and leading_dimension.  Otherwise
   // there is a problem with Mx1 and 1xM matrices.
   if (stride2(M) == 1) return 'T';
   if (stride1(M) == 1) return 'N';
   return 'X';
}

template <typename T>
inline
char blas_trans_row(T const& M)
{
   if (stride2(M) == 1) return 'N';
   if (stride1(M) == 1) return 'T';
   return 'X';
}

template <typename T>
inline
size_type leading_dimension(T const& M)
{
   if (stride2(M) == 1) return stride1(M);
   return stride2(M);
}


// is_row_major and is_col_major

template <typename T, typename Tinterface = typename interface<T>::type>
struct matrix_orientation;

template <typename T, typename Tv, typename Orient, typename Ti>
struct matrix_orientation<T, Concepts::DenseMatrix<Tv, Orient, Ti>>
{
   typedef Orient type;
};


template <typename T, typename Tinterface = typename interface<T>::type>
struct is_row_major : boost::mpl::false_ {};

template <typename T, typename Tv, typename Ti>
struct is_row_major<T, Concepts::DenseMatrix<Tv, RowMajor, Ti>> : boost::mpl::true_
{};

template <typename T, typename Tinterface = typename interface<T>::type>
struct is_col_major : boost::mpl::false_ {};

template <typename T, typename Tv, typename Ti>
struct is_col_major<T, Concepts::DenseMatrix<Tv, ColMajor, Ti>> : boost::mpl::true_
{};


} // namespace LinearAlgebra

#endif
