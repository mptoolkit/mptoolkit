/* -*- C++ -*- $Id$

  matrixsection.h

  proxy type for a matrix indirected by row & column projections.

  Created 2005-02-28 Ian McCulloch
*/

#if !defined(MATRIXSECTION_H_CHSJKH4UIWHUIHNLUIREHT8793Y8794YWFE80)
#define MATRIXSECTION_H_CHSJKH4UIWHUIHNLUIREHT8793Y8794YWFE80

#include "matrixoperations.h"
#include "matrixsectioniterator.h"
#include "slice.h"
#include "crtp_matrix.h"
#include <boost/mpl/assert.hpp>

namespace LinearAlgebra
{

template <typename MatType, typename RowIndex, typename ColIndex>
class MatrixSection : public MatrixBase<MatrixSection<MatType, RowIndex, ColIndex> >
{
   public:
      BOOST_MPL_ASSERT((is_proxy_reference<MatType>));
      BOOST_MPL_ASSERT((is_const_proxy_reference<RowIndex>));
      BOOST_MPL_ASSERT((is_const_proxy_reference<ColIndex>));

      typedef is_mutable_proxy_reference<MatType> proxy;
      typedef is_const_proxy_reference<MatType> const_proxy;

      typedef MatType matrix_reference;
      typedef RowIndex row_type;
      typedef ColIndex col_type;

      typedef typename interface<MatType>::value_type value_type;
      typedef typename make_const_reference<MatType>::type const_matrix_reference;

      MatrixSection(matrix_reference m, row_type r, col_type c)
         : m_(m), r_(r), c_(c) 
      {
         using LinearAlgebra::size1; using LinearAlgebra::size2;
	 DEBUG_PRECONDITION(size(r) == 0 || min(r) >= 0)(r);
	 DEBUG_PRECONDITION(size(r) == 0 || size_type(max(r)) < size1(m))(r)(size1(m));
	 DEBUG_PRECONDITION(size(c) == 0 || min(c) >= 0)(c); 
	 DEBUG_PRECONDITION(size(c) == 0 || size_type(max(c)) < size2(m))(c)(size2(m)); 
      }

      MatrixSection& operator=(MatrixSection const& x)
      {
         assign(*this, x);
         return *this;
      }

      template <typename U>
      typename boost::enable_if<is_matrix<U>, MatrixSection&>::type
      operator=(U const& x)
      {
         assign_copy(*this, x);
	 return *this;
      }

      template <typename U>
      typename boost::enable_if<is_matrix<U>, MatrixSection&>::type
      operator=(NoAliasProxy<U> const& x)
      {
	 assign(*this, x.value());
	 return *this;
      }

      size_type size1() const { return size(r_); }
      size_type size2() const { return size(c_); }

      const_matrix_reference base() const { return m_; }
      matrix_reference base() { return m_; }

      row_type row_index() const { return r_; }
      col_type col_index() const { return c_; }

   private:
      matrix_reference m_;
      row_type r_;
      col_type c_;
};

// interface

template <typename MatType, typename RowType, typename ColType,
          typename MatI = typename interface<MatType>::type,
          typename RowI = typename interface<RowType>::type,
          typename ColI = typename interface<ColType>::type>
struct MatrixSectionInterface {};

template <typename M, typename R, typename C,
          typename Mv, typename Orient, typename Mi,
          typename Rv, typename Ri,
          typename Cv, typename Ci>
struct MatrixSectionInterface<M, R, C,
                              DENSE_MATRIX(Mv, Orient, Mi),
                              DENSE_VECTOR(Rv, Ri),
                              DENSE_VECTOR(Cv, Ci)>
{
   typedef Mv value_type;
   typedef Orient orientation;
   typedef DENSE_MATRIX(Mv, Orient, void) type;
};

template <typename M,
          typename Mv, typename Orient, typename Mi>
struct MatrixSectionInterface<M, Range, Range,
                              STRIDE_MATRIX(Mv, Orient, Mi),
                              DENSE_VECTOR(size_type, void),
                              DENSE_VECTOR(size_type, void)>
{
   typedef Mv value_type;
   typedef Orient orientation;
   typedef STRIDE_MATRIX(Mv, Orient, void) type;
};

template <typename M,
          typename Mv, typename Orient, typename Mi>
struct MatrixSectionInterface<M, Range, Slice,
                              STRIDE_MATRIX(Mv, Orient, Mi),
                              DENSE_VECTOR(size_type, void),
                              DENSE_VECTOR(size_type, void)>
{
   typedef Mv value_type;
   typedef Orient orientation;
   typedef STRIDE_MATRIX(Mv, Orient, void) type;
};

template <typename M,
          typename Mv, typename Orient, typename Mi>
struct MatrixSectionInterface<M, Slice, Range,
                              STRIDE_MATRIX(Mv, Orient, Mi),
                              DENSE_VECTOR(size_type, void),
                              DENSE_VECTOR(size_type, void)>
{
   typedef Mv value_type;
   typedef Orient orientation;
   typedef STRIDE_MATRIX(Mv, Orient, void) type;
};

template <typename M,
          typename Mv, typename Orient, typename Mi>
struct MatrixSectionInterface<M, Slice, Slice,
                              STRIDE_MATRIX(Mv, Orient, Mi),
                              DENSE_VECTOR(size_type, void),
                              DENSE_VECTOR(size_type, void)>
{
   typedef Mv value_type;
   typedef Orient orientation;
   typedef STRIDE_MATRIX(Mv, Orient, void) type;
};

// FIXME: finish the interfaces, especially for sparse

template <typename M, typename R, typename C>
struct interface<MatrixSection<M, R, C> > : MatrixSectionInterface<M,R,C> {};

//
// TODO: the following should be specialized to only work if the MatrixSection
// actually conforms to the inferface.
//

template <typename M, typename RowIndex, typename ColIndex>
struct Stride1<MatrixSection<M, RowIndex, ColIndex> >
{
   typedef difference_type result_type;
   typedef MatrixSection<M, RowIndex, ColIndex> const& argument_type;

   result_type operator()(argument_type x) const
   {
      return stride1(x.base()) * stride(x.row_index());
   }
};

template <typename M, typename RowIndex, typename ColIndex>
struct Stride2<MatrixSection<M, RowIndex, ColIndex> >
{
   typedef difference_type result_type;
   typedef MatrixSection<M, RowIndex, ColIndex> const& argument_type;

   result_type operator()(argument_type x) const
   {
      return stride2(x.base()) * stride(x.col_index());
   }
};

template <typename M, typename RowIndex, typename ColIndex>
struct Data<MatrixSection<M, RowIndex, ColIndex> >
{
   typedef Data<typename basic_type<M>::type> Forward;
   typedef typename Forward::result_type result_type;
   typedef MatrixSection<M, RowIndex, ColIndex> const& argument_type;

   result_type operator()(argument_type x) const
   {
      return &x(0,0);
   }
};

template <typename M, typename RowIndex, typename ColIndex>
struct Data<MatrixSection<M, RowIndex, ColIndex>&>
{
   typedef Data<typename reference_to_arg<M>::type> Forward;
   typedef typename Forward::result_type result_type;
   typedef MatrixSection<M, RowIndex, ColIndex>& argument_type;

   result_type operator()(argument_type x) const
   {
      return &x(0,0);
   }
};

// iterators

template <typename M, typename R, typename C, 
          typename Orient = typename interface<MatrixSection<M,R,C> >::orientation>
struct MatrixSectionIterConst;

template <typename M, typename R, typename C>
struct MatrixSectionIterConst<M, R, C, RowMajor>
{
   typedef typename basic_type<M>::type MArg;
   typedef MatrixSectionOuterIterator<typename Iterate<MArg>::result_type,
                                      RowMajor,
                                      typename Iterate<typename basic_type<R>::type>::result_type,
                                      C>
      result_type;

   typedef MatrixSection<M, R, C> const& argument_type;
   result_type operator()(argument_type m) const
   {
      return result_type(iterate(m.base()), iterate(m.row_index()), m.col_index());
   }
};

template <typename M, typename R, typename C>
struct MatrixSectionIterConst<M, R, C, ColMajor>
{
   typedef typename basic_type<M>::type MArg;
   typedef MatrixSectionOuterIterator<typename Iterate<MArg>::result_type,
                                      ColMajor,
                                      typename Iterate<typename basic_type<C>::type>::result_type,
                                      R>
      result_type;

   typedef MatrixSection<M, R, C> const& argument_type;
   result_type operator()(argument_type m) const
   {
      return result_type(iterate(m.base()), iterate(m.col_index()), m.row_index());
   }
};

template <typename M, typename R, typename C, 
          typename Orient = typename interface<MatrixSection<M,R,C> >::orientation>
struct MatrixSectionIterMutable;

template <typename M, typename R, typename C>
struct MatrixSectionIterMutable<M, R, C, RowMajor>
{
   typedef typename reference_to_arg<M>::type MArg;
   //   typedef typename boost::mpl::print<M>::type dummy;
   //   typedef typename boost::mpl::print<MArg>::type dummy;
   //   typedef typename boost::mpl::print<typename Iterate<MArg>::result_type>::type dummy;
   typedef MatrixSectionOuterIterator<typename Iterate<MArg>::result_type,
                                      RowMajor,
                                      typename Iterate<typename basic_type<R>::type>::result_type,
                                      C>
      result_type;

   typedef MatrixSection<M, R, C>& argument_type;
   result_type operator()(argument_type m) const
   {
      return result_type(iterate(m.base()), iterate(m.row_index()), m.col_index());
   }
};

template <typename M, typename R, typename C>
struct MatrixSectionIterMutable<M, R, C, ColMajor>
{
   typedef typename reference_to_arg<M>::type MArg;
   typedef MatrixSectionOuterIterator<typename Iterate<MArg>::result_type,
                                      ColMajor,
                                      typename Iterate<typename basic_type<C>::type>::result_type,
                                      R>
      result_type;

   typedef MatrixSection<M, R, C>& argument_type;
   result_type operator()(argument_type m) const
   {
      return result_type(iterate(m.base()), iterate(m.col_index()), m.row_index());
   }
};

template <typename M, typename R, typename C>
struct Iterate<MatrixSection<M, R, C>&>
   : MatrixSectionIterMutable<M, R, C> {};

template <typename M, typename R, typename C>
struct Iterate<MatrixSection<M, R, C> >
   : MatrixSectionIterConst<M, R, C> {};

// size

// GetElement

template <typename M, typename R, typename C>
struct GetMatrixElement<MatrixSection<M, R, C>&,
                        typename boost::enable_if<
   boost::mpl::and_<
      is_mutable_proxy_reference<M>
    , is_defined<GetMatrixElement<typename reference_to_arg<M>::type> >
  >
>::type>
{
   typedef typename GetMatrixElement<typename 
   reference_to_arg<M>::type>::result_type result_type;

   typedef MatrixSection<M, R, C>& first_argument_type;
   typedef size_type second_argument_type;
   typedef size_type third_argument_type;

   result_type operator()(first_argument_type m, size_type i, size_type j) const
   {
      return get_element(m.base(), 
                         get_element(m.row_index(), i),
                         get_element(m.col_index(), j));
   }
};

template <typename M, typename R, typename C>
struct GetMatrixElement<MatrixSection<M, R, C>,
                        typename boost::enable_if<
   is_defined<GetMatrixElement<typename basic_type<M>::type> >
>::type>
{
   typedef typename GetMatrixElement<typename basic_type<M>::type>::result_type result_type;

   //typedef typename boost::mpl::print<typename basic_type<M>::type>::type dummy1;
   //typedef typename boost::mpl::print<M>::type dummy;

   typedef MatrixSection<M, R, C> const& first_argument_type;
   typedef size_type second_argument_type;
   typedef size_type third_argument_type;

   result_type operator()(first_argument_type m, size_type i, size_type j) const
   {
      return get_element(m.base(), 
                         get_element(m.row_index(), i),
                         get_element(m.col_index(), j));
   }
};

// SwapSortOrder

template <typename M, typename R, typename C>
struct SwapSortOrder<MatrixSection<M, R, C> >
{
   typedef typename SwapSortOrder<
      typename basic_type<M>::type
   >::result_type SwappedBase;

   BOOST_MPL_ASSERT((is_proxy_reference<SwappedBase>));

   typedef typename make_const_reference<SwappedBase>::type SwappedC;
   typedef MatrixSection<SwappedC, R, C> result_type;
   typedef MatrixSection<M, R, C> const& argument_type;

   result_type operator()(argument_type x) const 
   { return result_type(swap_sort_order(x.base()),x.row_index(), x.col_index()); }
};

template <typename M, typename R, typename C>
struct SwapSortOrder<MatrixSection<M, R, C>&>
{
   typedef typename SwapSortOrder<
      typename reference_to_arg<M>::type
   >::result_type SwappedBase;

   BOOST_MPL_ASSERT((is_proxy_reference<SwappedBase>));

   typedef MatrixSection<SwappedBase, R, C> result_type;
   typedef MatrixSection<M, R, C>& argument_type;

   result_type operator()(argument_type x) const 
   { return result_type(swap_sort_order(x.base()),x.row_index(), x.col_index()); }
};

//
// ProjectMatrix
//

template <typename T, typename U, typename V,
          typename Tv, typename Ti,
          typename Uv, typename Ui,
          typename Vv, typename Vi>
struct ProjectMatrixInterface<T, U, V,
                              ANY_MATRIX(Tv, Ti),
                              LOCAL_VECTOR(Uv, Ui),
                              LOCAL_VECTOR(Vv, Vi)>
{
   typedef MatrixSection<
      typename make_const_reference<T>::type,
      typename make_const_reference<U>::type,
      typename make_const_reference<V>::type
   > result_type;

   typedef T const& first_argument_type;
   typedef U const& second_argument_type;
   typedef V const& third_argument_type;

   result_type operator()(T const& m, U const& u, V const& v) const
   {
      return result_type(m,u,v);
   }
};

template <typename T, typename U, typename V,
          typename Tv, typename Ti,
          typename Uv, typename Ui,
          typename Vv, typename Vi>
struct ProjectMatrixInterface<T&, U, V,
                              ANY_MATRIX(Tv, Ti),
                              LOCAL_VECTOR(Uv, Ui),
                              LOCAL_VECTOR(Vv, Vi)>
{
   typedef MatrixSection<
      typename make_reference<T>::type,
      typename make_const_reference<U>::type,
      typename make_const_reference<V>::type
   > result_type;

   typedef T& first_argument_type;
   typedef U const& second_argument_type;
   typedef V const& third_argument_type;

   result_type operator()(T& m, U const& u, V const& v) const
   {
      return result_type(m,u,v);
   }
};

// MatrixRow,MatrixCol for stride matrices

template <typename T, typename Tv, typename Orient, typename Ti>
struct MatrixRowInterface<T, STRIDE_MATRIX(Tv, Orient, Ti)>
{
   typedef VectorMemProxy<Tv const, tagVariable> result_type;
   typedef T const& first_argument_type;
   typedef size_type second_argument_type;

   result_type operator()(first_argument_type v, second_argument_type n) const
   {
      DEBUG_RANGE_CHECK_OPEN(n, 0U, size1(v));
      return result_type(data(v) + n * stride1(v), size2(v), stride2(v));
   }
};

template <typename T, typename Tv, typename Orient, typename Ti>
struct MatrixRowInterface<T&, STRIDE_MATRIX(Tv, Orient, Ti)>
{
   typedef VectorMemProxy<Tv, tagVariable> result_type;
   typedef T& first_argument_type;
   typedef size_type second_argument_type;

   result_type operator()(first_argument_type v, second_argument_type n) const
   {
      DEBUG_RANGE_CHECK_OPEN(n, 0U, size1(v));
      return result_type(data(v) + n * stride1(v), size2(v), stride2(v));
   }
};

template <typename T, typename Tv, typename Orient, typename Ti>
struct MatrixColInterface<T, STRIDE_MATRIX(Tv, Orient, Ti)>
{
   typedef VectorMemProxy<Tv const, tagVariable> result_type;
   typedef T const& first_argument_type;
   typedef size_type second_argument_type;

   result_type operator()(first_argument_type v, second_argument_type n) const
   {
      DEBUG_RANGE_CHECK_OPEN(n, 0U, size2(v));
      return result_type(data(v) + n * stride2(v), size1(v), stride1(v));
   }
};

template <typename T, typename Tv, typename Orient, typename Ti>
struct MatrixColInterface<T&, STRIDE_MATRIX(Tv, Orient, Ti)>
{
   typedef VectorMemProxy<Tv, tagVariable> result_type;
   typedef T& first_argument_type;
   typedef size_type second_argument_type;

   result_type operator()(first_argument_type v, second_argument_type n) const
   {
      DEBUG_RANGE_CHECK_OPEN(n, 0U, size2(v));
      return result_type(data(v) + n * stride2(v), size1(v), stride1(v));
   }
};

} // namespace LinearAlgebra

#endif

