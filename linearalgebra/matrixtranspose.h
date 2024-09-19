// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/matrixtranspose.h
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
  matrixtranspose.h

  Created 2005-01-23 Ian McCulloch
*/

#if !defined(MATRIXTRANSPOSE_H_JDSHCEUIRHF735YT734YO87)
#define MATRIXTRANSPOSE_H_JDSHCEUIRHF735YT734YO87

#include "matrixiterators.h"
#include "matrixoperationsbase.h"
#include "crtp_matrix.h"

namespace LinearAlgebra
{

//
// iterator adaptors for the transpose
//

template <typename Iter>
class MatrixTransposedOuterIter : public IterAdaptorBase<Iter, MatrixTransposedOuterIter<Iter> >
{
   public:
      typedef IterAdaptorBase<Iter, MatrixTransposedOuterIter<Iter> > base_type;

#if 0
      typedef typename Iter::category category;
      typedef typename Iter::reference reference;
      typedef typename Iter::pointer pointer;
#endif

      MatrixTransposedOuterIter() {}

      MatrixTransposedOuterIter(Iter const& i) : base_type(i) {}

      size_type index() const { return base_type::index(); }

   using base_type::base;
   //      base_iterator const& base() const { return ; }
   //      base_iterator& base() { return I_; }


};

template <typename Iter>
class MatrixTransposedInnerIter : public IterAdaptorBase<Iter, MatrixTransposedInnerIter<Iter> >
{
   public:
      typedef IterAdaptorBase<Iter, MatrixTransposedInnerIter<Iter> > base_type;
      MatrixTransposedInnerIter() {}

      typedef typename Iter::category category;
      typedef typename Iter::reference reference;
      typedef typename Iter::pointer pointer;

      MatrixTransposedInnerIter(Iter const& i) : base_type(i) {}

      size_type index1() const { return this->base().index2(); }
      size_type index2() const { return this->base().index1(); }
};

template <typename Iter>
struct Iterate<MatrixTransposedOuterIter<Iter>&>
{
   typedef MatrixTransposedInnerIter<typename Iterate<Iter&>::result_type> result_type;
   typedef MatrixTransposedOuterIter<Iter> argument_type;
   result_type operator()(argument_type const& x) const
   {
      return result_type(iterate(x.base()));
   }
};

template <typename Iter>
struct Iterate<MatrixTransposedOuterIter<Iter> >
{
   typedef MatrixTransposedInnerIter<typename Iterate<Iter>::result_type> result_type;
   typedef MatrixTransposedOuterIter<Iter> argument_type;
   result_type operator()(argument_type const& x) const
   {
      return result_type(iterate(x.base()));
   }
};

template <typename BaseProxyReference>
class MatrixTransposeProxy : public MatrixBase<MatrixTransposeProxy<BaseProxyReference> >
{
   public:
      BOOST_MPL_ASSERT((is_proxy_reference<BaseProxyReference>));

      using MatrixBase<MatrixTransposeProxy<BaseProxyReference> >::operator();

      typedef is_const_proxy_reference<BaseProxyReference>   const_proxy;
      typedef is_mutable_proxy_reference<BaseProxyReference> proxy;

      typedef BaseProxyReference base_reference;
      typedef typename make_const_reference<BaseProxyReference>::type base_const_reference;
      typedef typename basic_type<BaseProxyReference>::type base_type;

      typedef typename Iterate<typename reference_to_arg<BaseProxyReference>::type>::result_type base_iterator;
      typedef typename const_iterator<base_type>::type const_base_iterator;

      typedef typename interface<base_type>::value_type value_type;
   //      typedef typename base_type::reference reference;
      typedef typename make_const_reference<base_reference>::type const_reference;

      typedef MatrixTransposedOuterIter<base_iterator> iterator;
      typedef MatrixTransposedOuterIter<const_base_iterator> const_iterator;

   //      typedef typename boost::remove_const<>::type mb_type;
   //      typedef typename boost::add_reference<mb_type>::type mutable_type;

      MatrixTransposeProxy(base_reference Base) : Base_(Base) {}

   // do we want transpose to be an l-value?
#if 0
      reference operator()(size_type i, size_type j)
         { return Base_(j,i); }
#endif

      const_reference operator()(size_type i, size_type j) const
         { return Base_(j,i); }

      base_reference base() { return Base_; }

      base_const_reference base() const { return Base_; }
      base_reference mutable_base() const { return Base_; }

   private:
      base_reference Base_;
};

// basic operations

template <typename M>
struct Size1<MatrixTransposeProxy<M> >
{
   typedef typename basic_type<M>::type base_type;
   typedef typename Size2<base_type>::result_type result_type;
   typedef MatrixTransposeProxy<M> const& argument_type;
   result_type operator()(argument_type m) const
   {
      return Size2<base_type>()(m.base());
   }
};

template <typename M>
struct Size2<MatrixTransposeProxy<M> >
{
   typedef typename basic_type<M>::type base_type;
   typedef typename Size1<base_type>::result_type result_type;
   typedef MatrixTransposeProxy<M> const& argument_type;
   result_type operator()(argument_type m) const
   {
      return Size1<base_type>()(m.base());
   }
};

template <typename M>
struct Stride1<MatrixTransposeProxy<M> >
{
   typedef typename basic_type<M>::type base_type;
   typedef typename Stride2<base_type>::result_type result_type;
   typedef MatrixTransposeProxy<M> const& argument_type;
   result_type operator()(argument_type m) const
   {
      return Stride2<base_type>()(m.base());
   }
};

template <typename M>
struct Stride2<MatrixTransposeProxy<M> >
{
   typedef typename basic_type<M>::type base_type;
   typedef typename Stride1<base_type>::result_type result_type;
   typedef MatrixTransposeProxy<M> const& argument_type;
   result_type operator()(argument_type m) const
   {
      return Stride1<base_type>()(m.base());
   }
};
template <typename M>
struct Data<MatrixTransposeProxy<M> >
{
   typedef Data<typename basic_type<M>::type> Fwd;
   typedef typename Fwd::result_type result_type;
   typedef MatrixTransposeProxy<M> argument_type;
   result_type operator()(MatrixTransposeProxy<M> const& x) const
   { return data(x.base()); }
};

template <typename M>
struct Data<MatrixTransposeProxy<M>&>
{
   typedef Data<typename reference_to_arg<M>::type> Fwd;
   typedef typename Fwd::result_type result_type;
   typedef MatrixTransposeProxy<M>& argument_type;
   result_type operator()(MatrixTransposeProxy<M>& x) const
   { return data(x.base()); }
};

// interface

template <typename T, typename NewT, typename Ti = typename interface<T>::type>
struct TransposeOrientationInterface : interface<T> {};

template <typename T, typename NewT, typename Tv, typename Orient, typename Ti>
struct TransposeOrientationInterface<T, NewT, Concepts::CompressedOuterMatrix<Tv, Orient, Ti>>
{
   typedef typename SwapOrientation<Orient>::type NewOrient;
   typedef Concepts::CompressedOuterMatrix<Tv, NewOrient, NewT> type;
   typedef Tv value_type;
};

template <typename T, typename NewT, typename Tv, typename Orient, typename Ti>
struct TransposeOrientationInterface<T, NewT, Concepts::DenseMatrix<Tv, Orient, Ti>>
{
   typedef typename SwapOrientation<Orient>::type NewOrient;
   typedef Concepts::DenseMatrix<Tv, NewOrient, NewT> type;
   typedef Tv value_type;
};

template <typename T, typename NewT, typename Tv, typename Orient, typename Ti>
struct TransposeOrientationInterface<T, NewT, Concepts::StrideMatrix<Tv, Orient, Ti>>
{
   typedef typename SwapOrientation<Orient>::type NewOrient;
   typedef Concepts::StrideMatrix<Tv, NewOrient, NewT> type;
   typedef Tv value_type;
};

template <typename T, typename NewT, typename Tv, typename Orient, typename Ti>
struct TransposeOrientationInterface<T, NewT, Concepts::ContiguousMatrix<Tv, Orient, Ti>>
{
   typedef typename SwapOrientation<Orient>::type NewOrient;
   typedef Concepts::ContiguousMatrix<Tv, NewOrient, NewT> type;
   typedef Tv value_type;
};


template <typename M>
struct interface<MatrixTransposeProxy<M> >
   : TransposeOrientationInterface<M, void > {};

// iterators

template <typename M>
struct Iterate<MatrixTransposeProxy<M>&>
{
   typedef typename MatrixTransposeProxy<M>::iterator result_type;
   typedef MatrixTransposeProxy<M>& argument_type;
   result_type operator()(argument_type x) const
   { return result_type(iterate(x.base())); }
};

template <typename M>
struct Iterate<MatrixTransposeProxy<M> >
{
   typedef typename MatrixTransposeProxy<M>::const_iterator result_type;
   typedef MatrixTransposeProxy<M> const& argument_type;
   result_type operator()(argument_type x) const
   { return result_type(iterate(x.base())); }
};

// row/col access

template <typename M, typename Enable = void>
struct MatrixRow_MatrixTransposeProxy {};

template <typename M>
struct MatrixRow_MatrixTransposeProxy<
   M,
   typename boost::enable_if<exists<typename MatrixCol<M>::result_type> >::type>
{
   typedef typename MatrixCol<M>::result_type result_type;
   typedef MatrixTransposeProxy<M> first_argument_type;
   typedef size_type second_argument_type;
   result_type operator()(first_argument_type x,
      size_type n) const { return matrix_col(x, n); }
};

template <typename M>
struct MatrixRow<MatrixTransposeProxy<M> > : MatrixRow_MatrixTransposeProxy<M> {};

template <typename M, typename Enable = void>
struct MatrixCol_MatrixTransposeProxy {};

template <typename M>
struct MatrixCol_MatrixTransposeProxy<
   M,
   typename boost::enable_if<exists<typename MatrixRow<M>::result_type> >::type>
{
   typedef typename MatrixRow<M>::result_type result_type;
   typedef MatrixTransposeProxy<M> first_argument_type;
   typedef size_type second_argument_type;
   result_type operator()(first_argument_type x,
      size_type n) const { return matrix_row(x, n); }
};

template <typename M>
struct MatrixCol<MatrixTransposeProxy<M> > : MatrixCol_MatrixTransposeProxy<M> {};

// transpose

template <typename M, typename F, typename Enable>
struct MatrixTranspose<MatrixTransposeProxy<M>, F, Enable>
{
   typedef boost::mpl::true_ involutary;
   typedef typename Transform<M, F>::result_type result_type;
   typedef MatrixTransposeProxy<M> const& argument_type;
   typedef argument_type first_argument_type;
   typedef F second_argument_type;
   result_type operator()(argument_type x) const { return transform(x.base(), F()); }
   result_type operator()(first_argument_type x,
                          second_argument_type f) const { return transform(x.base(), f); }
};

template <typename M, typename F, typename Enable>
struct MatrixTranspose<MatrixTransposeProxy<M>&, F, Enable>
{
   typedef boost::mpl::true_ involutary;
   typedef typename Transform<M&, F>::result_type result_type;
   typedef MatrixTransposeProxy<M>& argument_type;
   typedef argument_type first_argument_type;
   typedef F second_argument_type;
   result_type operator()(argument_type x) const { return transform(x.base(), F()); }
   result_type operator()(first_argument_type x,
                          second_argument_type f) const { return transform(x.base(), f); }
};



template <typename M, typename F>
struct MatrixTranspose<MatrixTransposeProxy<M>, F, typename boost::enable_if<is_identity<F> >::type>
{
   typedef boost::mpl::true_ involutary;
   typedef M const result_type;
   typedef MatrixTransposeProxy<M> const& argument_type;
   typedef argument_type first_argument_type;
   typedef F second_argument_type;
   result_type operator()(argument_type x) const { return x.base(); }
   result_type operator()(first_argument_type x,
                          second_argument_type const&) const { return x.base(); }
};

template <typename M, typename F>
struct MatrixTranspose<MatrixTransposeProxy<M>&, F, typename boost::enable_if<is_identity<F> >::type>
{
   typedef boost::mpl::true_ involutary;
   typedef M result_type;
   typedef MatrixTransposeProxy<M>& argument_type;
   typedef argument_type first_argument_type;
   typedef F second_argument_type;
   result_type operator()(argument_type x) const { return x.base(); }
   result_type operator()(first_argument_type x,
                          second_argument_type const&) { return x.base(); }
};

template <typename M, typename F, typename Mv, typename Mi, typename Enable>
struct MatrixTransposeInterface<M, F, Concepts::LocalMatrix<Mv, Mi>, Enable>
{
   typedef boost::mpl::true_ involutary;
   typedef MatrixTransposeProxy<typename Transform<M, F>::result_type> result_type;
   typedef M const& argument_type;
   typedef argument_type first_argument_type;
   typedef F second_argument_type;
   result_type operator()(M const& m) const { return result_type(transform(m, F())); }
   result_type operator()(first_argument_type m,
                          second_argument_type f) const { return result_type(transform(m,f)); }
};

template <typename M, typename F, typename Mv, typename Mi>
struct MatrixTransposeInterface<M&, F, Concepts::LocalMatrix<Mv, Mi>>
{
   typedef boost::mpl::true_ involutary;
   typedef MatrixTransposeProxy<typename Transform<M&, F>::result_type> result_type;
   typedef M& argument_type;
   typedef argument_type first_argument_type;
   typedef F second_argument_type;
   result_type operator()(M& m) const { return result_type(transform(m, F())); }
   result_type operator()(first_argument_type m,
                          second_argument_type f) const { return result_type(transform(m,f)); }
};

// SwapSortOrder

template <typename M>
struct SwapSortOrder<MatrixTransposeProxy<M> >
{
   typedef boost::mpl::true_ involutary;
   typedef typename reference_to_arg<M>::type MArg;
   typedef typename SwapSortOrder<MArg>::result_type SwappedBase;
   typedef typename Transpose<typename reference_to_arg<SwappedBase>::type>::result_type
   result_type;
   typedef MatrixTransposeProxy<M> const& argument_type;
   result_type operator()(argument_type x) const
   { return transpose(swap_sort_order(x.base())); }
};

template <typename M>
struct SwapSortOrder<MatrixTransposeProxy<M>&>
{
   typedef boost::mpl::true_ involutary;
   typedef typename reference_to_arg<M>::type MArg;
   typedef typename SwapSortOrder<MArg>::result_type SwappedBase;
   typedef typename Transpose<typename reference_to_arg<SwappedBase>::type>::result_type result_type;
   typedef MatrixTransposeProxy<M>& argument_type;
   result_type operator()(argument_type x) const
   { return transpose(swap_sort_order(x.base())); }
};

// expression operations

} // namespace LinearAlgebra

#endif
