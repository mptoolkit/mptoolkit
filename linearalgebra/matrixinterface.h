/* -*- C++ -*- $Id$

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

template <typename Value, typename T>
struct MatrixExpression
{
};

#define ANY_MATRIX(Value, T) LinearAlgebra::Concepts::MatrixExpression<Value, T > /**/
#define MATRIX_EXPRESSION(Value, T) LinearAlgebra::Concepts::MatrixExpression<Value, T > /**/

#define LOCAL_MATRIX(Value, T)				\
LinearAlgebra::Concepts::MatrixExpression< Value ,	\
   LinearAlgebra::Concepts::Local< T > > /**/


template <typename T>
struct SparseMatrix {}; // : Local<T> {};

#define SPARSE_MATRIX(Value, T)					\
LinearAlgebra::Concepts::MatrixExpression< Value ,		\
   LinearAlgebra::Concepts::Local<				\
      LinearAlgebra::Concepts::SparseMatrix< T > > > /**/

template <typename T>
struct Rectangular {}; // : SparseMatrix<Rectangular<T> > {};

#define RECTANGULAR_MATRIX(Value, T)				\
LinearAlgebra::Concepts::MatrixExpression< Value ,		\
   LinearAlgebra::Concepts::Local<				\
      LinearAlgebra::Concepts::SparseMatrix<			\
         LinearAlgebra::Concepts::Rectangular< T > > > /**/


#define COMPRESSED_MATRIX(Value, T)				\
LinearAlgebra::Concepts::MatrixExpression< Value ,		\
   LinearAlgebra::Concepts::Local<				\
      LinearAlgebra::Concepts::SparseMatrix<			\
         LinearAlgebra::Concepts::Rectangular<			\
            LinearAlgebra::Concepts::Compressed< T > > > > > /**/


template <typename T>
struct Coordinate {}; // : Rectangular<Coordinate<T> > {};

#define COORDINATE_MATRIX(Value, T)				\
LinearAlgebra::Concepts::MatrixExpression< Value ,		\
   LinearAlgebra::Concepts::Local<				\
      LinearAlgebra::Concepts::SparseMatrix<			\
         LinearAlgebra::Concepts::Rectangular<			\
            LinearAlgebra::Concepts::Coordinate< T > > > > > /**/


template <typename T>
struct InjectiveCoordinate {}; // : Coordinate<InjectiveCoordinate<T> > {};

#define INJECTIVE_MATRIX(Value, T)						\
LinearAlgebra::Concepts::MatrixExpression< Value ,				\
   LinearAlgebra::Concepts::Local<						\
      LinearAlgebra::Concepts::SparseMatrix<					\
         LinearAlgebra::Concepts::Rectangular<					\
            LinearAlgebra::Concepts::Coordinate<				\
               LinearAlgebra::Concepts::InjectiveCoordinate< T > > > > > > /**/


template <typename T>
struct Banded {}; // : SparseMatrix<Banded<T> > {};

#define BANDED_MATRIX(Value, T)					\
LinearAlgebra::Concepts::MatrixExpression< Value ,		\
   LinearAlgebra::Concepts::Local<				\
      LinearAlgebra::Concepts::SparseMatrix<			\
         LinearAlgebra::Concepts::Banded< T > > > > /**/


template <typename T>
struct Diagonal {}; // : Banded<Diagonal<T> > {};

#define DIAGONAL_MATRIX(Value, T)				\
LinearAlgebra::Concepts::MatrixExpression< Value ,		\
   LinearAlgebra::Concepts::Local<				\
      LinearAlgebra::Concepts::SparseMatrix<			\
         LinearAlgebra::Concepts::Banded<			\
            LinearAlgebra::Concepts::Diagonal< T > > > > > /**/


template <typename Orientation, typename T>
struct CompressedOuter {}; // : Rectangular<CompressedOuter<Orientation, T> > {};

#define COMPRESSED_OUTER_MATRIX(Value, Orient, T)					\
LinearAlgebra::Concepts::MatrixExpression< Value ,					\
   LinearAlgebra::Concepts::Local<							\
      LinearAlgebra::Concepts::SparseMatrix<						\
         LinearAlgebra::Concepts::Rectangular<						\
            LinearAlgebra::Concepts::Compressed<					\
               LinearAlgebra::Concepts::CompressedOuter< Orient, T > > > > > > /**/

#define COMPRESSED_ROW_MATRIX(Value, T)							\
LinearAlgebra::Concepts::MatrixExpression< Value ,					\
   LinearAlgebra::Concepts::Local<							\
      LinearAlgebra::Concepts::SparseMatrix<						\
         LinearAlgebra::Concepts::Rectangular<						\
            LinearAlgebra::Concepts::Compressed<					\
            LinearAlgebra::Concepts::CompressedOuter< RowMajor, T > > > > > > /**/

#define COMPRESSED_COL_MATRIX(Value, T)							\
LinearAlgebra::Concepts::MatrixExpression< Value ,					\
   LinearAlgebra::Concepts::Local<							\
      LinearAlgebra::Concepts::SparseMatrix<						\
         LinearAlgebra::Concepts::Rectangular<						\
            LinearAlgebra::Concepts::Compressed<					\
               LinearAlgebra::Concepts::CompressedOuter< ColMajor, T > > > > > > /**/


template <typename Orientation, typename OuterInterface, typename T>
struct CompressedMatrix {};
   //   : CompressedOuter<Orientation, CompressedMatrix<Orientation, OuterInterface, T> > {};

#define COMPRESSED_OUTER_MATRIX_V(Value, Orient, VecInterface, T)	\
LinearAlgebra::Concepts::MatrixExpression< Value ,			\
   LinearAlgebra::Concepts::Local<					\
      LinearAlgebra::Concepts::SparseMatrix<				\
         LinearAlgebra::Concepts::Rectangular<				\
            LinearAlgebra::Concepts::Compressed<			\
               LinearAlgebra::Concepts::CompressedOuter< Orient,	\
                  LinearAlgebra::Concepts::CompressedMatrix< 		\
    Orient, VecInterface, T > > > > > > > /**/



template <typename Orientation, typename T>
struct DenseMatrix {}; // : Local<DenseMatrix<Orientation, T> > {};

#define DENSE_MATRIX(Value, Orient, T)				\
LinearAlgebra::Concepts::MatrixExpression< Value ,		\
   LinearAlgebra::Concepts::Local<				\
      LinearAlgebra::Concepts::DenseMatrix<Orient, T > > > /**/

#define DENSE_ROW_MAJOR_MATRIX(Value, T)				\
LinearAlgebra::Concepts::MatrixExpression< Value ,			\
   LinearAlgebra::Concepts::Local<					\
      LinearAlgebra::Concepts::DenseMatrix<RowMajor, T > > > /**/

#define DENSE_COL_MAJOR_MATRIX(Value, T)				\
LinearAlgebra::Concepts::MatrixExpression< Value ,			\
   LinearAlgebra::Concepts::Local<					\
      LinearAlgebra::Concepts::DenseMatrix<ColMajor, T > > > /**/


template <typename Orientation, typename T>
struct StrideMatrix {}; // : DenseMatrix<Orientation, StrideMatrix<Orientation, T> > {};

#define STRIDE_MATRIX(Value, Orient, T)						\
LinearAlgebra::Concepts::MatrixExpression< Value ,				\
   LinearAlgebra::Concepts::Local<						\
      LinearAlgebra::Concepts::DenseMatrix<Orient,				\
         LinearAlgebra::Concepts::StrideMatrix<Orient, T > > > > /**/

#define STRIDE_ROW_MAJOR_MATRIX(Value, T)					\
LinearAlgebra::Concepts::MatrixExpression< Value ,				\
   LinearAlgebra::Concepts::Local<						\
      LinearAlgebra::Concepts::DenseMatrix<RowMajor,				\
         LinearAlgebra::Concepts::StrideMatrix<RowMajor, T > > > > /**/

#define STRIDE_COL_MAJOR_MATRIX(Value, T)					\
LinearAlgebra::Concepts::MatrixExpression< Value ,				\
   LinearAlgebra::Concepts::Local<						\
      LinearAlgebra::Concepts::DenseMatrix<ColMajor,				\
         LinearAlgebra::Concepts::StrideMatrix<ColMajor, T > > > > /**/


template <typename Orientation, typename T>
struct ContiguousMatrix {}; // : StrideMatrix<Orientation, ContiguousMatrix<Orientation, T> > {};

#define CONTIGUOUS_MATRIX(Value, Orient, T)						\
LinearAlgebra::Concepts::MatrixExpression< Value ,					\
   LinearAlgebra::Concepts::Local<							\
      LinearAlgebra::Concepts::DenseMatrix<Orient,					\
         LinearAlgebra::Concepts::StrideMatrix<Orient,				\
            LinearAlgebra::Concepts::ContiguousMatrix<Orient, T > > > > > /**/

#define CONTIGUOUS_ROW_MAJOR_MATRIX(Value, T)						\
LinearAlgebra::Concepts::MatrixExpression< Value ,					\
   LinearAlgebra::Concepts::Local<							\
      LinearAlgebra::Concepts::DenseMatrix<RowMajor,					\
         LinearAlgebra::Concepts::StrideMatrix<RowMajor,				\
            LinearAlgebra::Concepts::ContiguousMatrix<RowMajor, T > > > > > /**/

#define CONTIGUOUS_COL_MAJOR_MATRIX(Value, T)						\
LinearAlgebra::Concepts::MatrixExpression< Value ,					\
   LinearAlgebra::Concepts::Local<							\
      LinearAlgebra::Concepts::DenseMatrix<ColMajor,					\
         LinearAlgebra::Concepts::StrideMatrix<ColMajor,				\
            LinearAlgebra::Concepts::ContiguousMatrix<ColMajor, T > > > > > /**/


} // namespace Concepts

// RebindInterface

template <typename T, typename Vv, typename Vi, typename Tv>
struct RebindInterface<T, MATRIX_EXPRESSION(Vv, Vi), Tv>
{
   typedef Tv value_type;
   typedef MATRIX_EXPRESSION(Tv, void) type;
};
template <typename T, typename Vv, typename Vi, typename Tv>
struct RebindInterface<T, LOCAL_MATRIX(Vv, Vi), Tv>
{
   typedef Tv value_type;
   typedef LOCAL_MATRIX(Tv, void) type;
};
template <typename T, typename Vv, typename Vi, typename Tv>
struct RebindInterface<T, COMPRESSED_MATRIX(Vv, Vi), Tv>
{
   typedef Tv value_type;
   typedef COMPRESSED_MATRIX(Tv, void) type;
};
template <typename T, typename Vv, typename Vi, typename Tv>
struct RebindInterface<T, COORDINATE_MATRIX(Vv, Vi), Tv>
{
   typedef Tv value_type;
   typedef COORDINATE_MATRIX(Tv, void) type;
};
template <typename T, typename Vv, typename Vi, typename Tv>
struct RebindInterface<T, INJECTIVE_MATRIX(Vv, Vi), Tv>
{
   typedef Tv value_type;
   typedef INJECTIVE_MATRIX(Tv, void) type;
};
template <typename T, typename Vv, typename Vi, typename Tv>
struct RebindInterface<T, BANDED_MATRIX(Vv, Vi), Tv>
{
   typedef Tv value_type;
   typedef BANDED_MATRIX(Tv, void) type;
};
template <typename T, typename Vv, typename Vi, typename Tv>
struct RebindInterface<T, DIAGONAL_MATRIX(Vv, Vi), Tv>
{
   typedef Tv value_type;
   typedef DIAGONAL_MATRIX(Tv, void) type;
};
template <typename T, typename Vv, typename Vo, typename Vi, typename Tv>
struct RebindInterface<T, COMPRESSED_OUTER_MATRIX(Vv, Vo, Vi), Tv>
{
   typedef Tv value_type;
   typedef COMPRESSED_OUTER_MATRIX(Tv, Vo, void) type;
};

#if 0
template <typename T, typename Vv, typename Vo, typename Vi, typename Tv>
struct RebindInterface<T, INJECTIVE_OUTER_MATRIX(Vv, Vo, Vi), Tv>
{
   typedef INJECTIVE_OUTER_MATRIX(Tv, Vo, void) type;
};
template <typename T, typename Vv, typename Vo, typename Vi, typename Tv>
struct RebindInterface<T, DENSE_OUTER_MATRIX(Vv, Vo, Vi), Tv>
{
   typedef DENSE_OUTER_MATRIX(Tv, Vo, void) type;
};
#endif

template <typename T, typename Vv, typename Vo, typename Vi, typename Tv>
struct RebindInterface<T, DENSE_MATRIX(Vv, Vo, Vi), Tv>
{
   typedef DENSE_MATRIX(Tv, Vo, void) type;
};
template <typename T, typename Vv, typename Vo, typename Vi, typename Tv>
struct RebindInterface<T, STRIDE_MATRIX(Vv, Vo, Vi), Tv>
{
   typedef STRIDE_MATRIX(Tv, Vo, void) type;
};
template <typename T, typename Vv, typename Vo, typename Vi, typename Tv>
struct RebindInterface<T, CONTIGUOUS_MATRIX(Vv, Vo, Vi), Tv>
{
   typedef CONTIGUOUS_MATRIX(Tv, Vo, void) type;
};

// is_matrix, boolean function to determine if a type
// has a matrix interface.

namespace Private
{
template <typename T>
struct is_matrix_helper : boost::mpl::false_ {};

template <typename S, typename T>
struct is_matrix_helper<MATRIX_EXPRESSION(S,T)> : boost::mpl::true_ {};

template <typename T>
struct is_dense_matrix_helper : boost::mpl::false_ {};

template <typename S, typename T, typename U>
struct is_dense_matrix_helper<DENSE_MATRIX(S,T,U)> : boost::mpl::true_ {};

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


//
// make_vector_from_abstract
//
// provides a 'default' choice of a concrete vector type based on
// the abstract type (ie. either sparse or dense).
//

template <typename T, typename Abstract>
struct make_matrix_from_abstract {};

template <typename T, typename Sv, typename So, typename Si>
struct abstract_interface_interface<T, DENSE_MATRIX(Sv, So, Si)>
   : matrix_abstract_dense {};

template <typename T, typename Sv, typename Si>
struct abstract_interface_interface<T, SPARSE_MATRIX(Sv, Si)>
   : matrix_abstract_sparse {};

// make_value_from_interface for vectors

template <typename T, typename Sv, typename Si>
struct make_value_from_interface<T, ANY_MATRIX(Sv, Si)>
   : make_matrix_from_abstract<Sv, typename abstract_interface<T>::type> {};


//
// iterator categories - unused?
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
struct matrix_orientation<T, DENSE_MATRIX(Tv, Orient, Ti)>
{
   typedef Orient type;
};


template <typename T, typename Tinterface = typename interface<T>::type>
struct is_row_major : boost::mpl::false_ {};

template <typename T, typename Tv, typename Ti>
struct is_row_major<T, DENSE_MATRIX(Tv, RowMajor, Ti)> : boost::mpl::true_ 
{};

template <typename T, typename Tinterface = typename interface<T>::type>
struct is_col_major : boost::mpl::false_ {};

template <typename T, typename Tv, typename Ti>
struct is_col_major<T, DENSE_MATRIX(Tv, ColMajor, Ti)> : boost::mpl::true_ 
{};


} // namespace LinearAlgebra

#endif
