// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/matrixtransform.h
//
// Copyright (C) 2005-2016 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2012 Stefan Depenbrock <Stefan.Depenbrock@physik.uni-muenchen.de>
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
  matrixtransform.h

  Created 2005-01-14 Ian McCulloch
*/

#if !defined(MATRIXTRANSFORM_H_DJSHC84YT789YUFP89RHJT89PJHFP8943)
#define MATRIXTRANSFORM_H_DJSHC84YT789YUFP89RHJT89PJHFP8943

#include "matrixinterface.h"
#include "matrixoperationsbase.h"
#include "matrixtranspose.h"
#include "matrixtransformiterator.h"
#include <boost/mpl/assert.hpp>

#include "matrix.h" // for a temp-matrix type

namespace LinearAlgebra
{

template <typename BaseProxyReference, typename F>
class MatrixTransformProxy
{
   public:

      typedef typename std::decay<typename std::remove_reference<BaseProxyReference>::type>::type BaseType;
      typedef typename std::result_of<F(typename BaseType::value_type)>::type reference;
      //      typedef typename F::result_type reference;
      typedef typename make_const_reference<reference>::type const_reference;
      typedef typename make_value<const_reference>::type value_type;

      typedef is_mutable_proxy_reference<reference> proxy;
      typedef boost::mpl::not_<proxy> const_proxy;

   //typedef typename boost::mpl::print<typename const_proxy::type>::type dummy;

      typedef BaseProxyReference base_reference;
   //      typedef typename make_const_reference<base_reference>::type base_const_reference;

       // note that we cannot define
       // typedef MatrixTransformProxy<base_const_reference, F> const_type;
       // here, because F may be expecting to take its argument type by non-const reference,
       // so we must use base_reference rather than base_const_reference.
       // Instead, we must use
       // typedef MatrixTransformProxy const const_type;
       // but this is the default, so we omit it.

      typedef F functor_type;

      // declare the abstract interface; only used if BaseType is an expression
      //typedef typename abstract_interface<BaseProxyReference>::type abstract_interface;

      MatrixTransformProxy(BaseProxyReference Base, functor_type const& Func)
        : Base_(Base), Func_(Func) {}

      explicit MatrixTransformProxy(BaseProxyReference Base)
        : Base_(Base), Func_() {}

      MatrixTransformProxy(MatrixTransformProxy const& m)
         : Base_(m.Base_), Func_(m.Func_) {}

      // this conversion ctor handles non-const to const conversions.
      template <typename OtherBase>
      MatrixTransformProxy(MatrixTransformProxy<OtherBase, F> const& e)
         : Base_(e.base()), Func_(e.functor()) {}

      size_type size1() const { using LinearAlgebra::size1; return size1(Base_); }
      size_type size2() const { using LinearAlgebra::size2; return size2(Base_); }

      const_reference operator()(size_type i, size_type j) const
         { return Func_(Base_(i,j)); }

      reference operator()(size_type i, size_type j)
         { return Func_(Base_(i,j)); }

      template <typename U>
      typename boost::enable_if<is_matrix<U>, MatrixTransformProxy&>::type
      operator=(U const& x)
      {
         // TODO: better temp type
         Matrix<value_type> Temp(x);

         assign(*this, Temp);
         return *this;
      }

      template <typename U>
      typename boost::enable_if<is_matrix<U>, MatrixTransformProxy&>::type
      operator=(NoAliasProxy<U> const& x)
      {
         assign(*this, x.value());
         return *this;
      }

      base_reference base() { return Base_; }
      base_reference base() const { return Base_; }

      functor_type const& functor() const { return Func_; }

      base_reference mutable_base() const { return Base_; }

   private:
      base_reference Base_;
      functor_type Func_;
};

// interface

template <typename T, typename OldInterface>
struct MatrixTransformInterface;

// generic version
template <typename T, typename U, typename V>
struct MatrixTransformInterface<T, Concepts::AnyMatrix<U,V>>
{
   // FIXME: this gets the Value type wrong, but we don't have RebindInterface anymore -
   // can we resolve that?  Do we really need the Value type to be part of the concept?
   typedef Concepts::AnyMatrix<U,V> type;
   //typedef typename RebindInterface<T, Concepts::AnyMatrix<U,V>>::type type;
};

// cut-off the interface rebinding at DENSE_MATRIX - we cannot have
// a MatrixTransform that is a stride or contiguous matrix.
template <typename T, typename U, typename Orient, typename V>
struct MatrixTransformInterface<T, Concepts::DenseMatrix<U,Orient,V>>
{
   typedef Concepts::DenseMatrix<typename T::value_type, Orient, void> type;
};

template <typename BaseType, typename F>
struct interface<MatrixTransformProxy<BaseType, F> >
{
   typedef typename MatrixTransformProxy<BaseType, F>::value_type value_type;
   typedef typename MatrixTransformInterface<
      MatrixTransformProxy<BaseType, F>, typename interface<BaseType>::type>::type type;
};

// for a MatrixTransformProxy of a MATRIX_EXPRESSION, we need to
// supply EvalExpression
template <typename Base, typename F, typename Value>
struct EvalExpression<MatrixTransformProxy<Base, F>,
                      Concepts::MatrixExpression<Value, MatrixTransformProxy<Base, F>>>
{
   typedef typename EvalExpression<Base>::result_type BaseValueType;
   typedef Transform<BaseValueType, F> Transformer;

   //   typedef typename Transformer::result_type result_type;
   typedef typename make_value<typename Transformer::result_type>::type result_type;
   typedef MatrixTransformProxy<Base, F> argument_type;

   result_type operator()(argument_type const& x) const
   {
      return result_type(transform(eval_expression(x.base()), x.functor()));
   }
};

#if 0
template <typename Base, typename F, typename Value>
struct EvalExpression<MatrixTransformProxy<Base, F> const,
                      Concepts::MatrixExpression<Value, MatrixTransformProxy<Base, F> const>>>
{
   typedef typename EvalExpression<Base>::result_type BaseValueType;
   typedef Transform<BaseValueType, F> Transformer;

   typedef typename make_value<typename Transformer::result_type>::type result_type;
   typedef MatrixTransformProxy<Base, F> argument_type;

   result_type operator()(argument_type const& x) const
   {
      return result_type(transform(eval_expression(x.base()), x.functor()));
   }
};
#endif

// iterators

template <typename Base, typename F>
struct Iterate<MatrixTransformProxy<Base, F>&>
{
   typedef MatrixTransformProxy<Base, F>& argument_type;
   typedef typename iterator<typename reference_to_arg<Base>::type>::type base_iterator;
   typedef MatrixTransformOuterIterator<base_iterator, F> result_type;
   result_type operator()(argument_type x) const
   {
      return result_type(iterate(x.base()), x.functor());
   }
};

template <typename Base, typename F>
struct Iterate<MatrixTransformProxy<Base, F> >
{
   typedef MatrixTransformProxy<Base, F> const& argument_type;
   typedef typename Iterate<typename reference_to_arg<Base>::type>::result_type base_iterator;

   //   typedef typename boost::mpl::print<base_iterator>::type dummy;

   typedef MatrixTransformOuterIterator<base_iterator, F> result_type;
   result_type operator()(argument_type x) const
   {
      return result_type(iterate(x.mutable_base()), x.functor());
   }
};

// transform

template <typename T, typename F>
struct TransformMatrix
{
   typedef T const& argument_type;
   typedef T const& first_argument_type;
   typedef F second_argument_type;

   typedef MatrixTransformProxy<typename make_const_reference<T>::type, F> result_type;

   result_type operator()(argument_type E) const
      { return result_type(E, F()); }

   result_type operator()(first_argument_type E, second_argument_type const& f) const
      { return result_type(E, f); }
};

template <typename T, typename F>
struct TransformMatrix<T&, F>
{
   typedef T& argument_type;
   typedef T& first_argument_type;
   typedef F const& second_argument_type;

   typedef MatrixTransformProxy<typename make_reference<T>::type, F> result_type;

   result_type operator()(argument_type E) const
      { return result_type(E, F()); }

   result_type operator()(first_argument_type E, second_argument_type f) const
      { return result_type(E, f); }
};

template <typename Base, typename F, typename G>
struct TransformMatrix<MatrixTransformProxy<Base, F>&, G>
{
   typedef MatrixTransformProxy<Base, F>& argument_type;
   typedef MatrixTransformProxy<Base, F>& first_argument_type;
   typedef G second_argument_type;
   typedef Compose<G, F> Composer;
   typedef typename Composer::result_type Func;

   typedef Transform<typename reference_to_arg<Base>::type, Func> Transformer;
   typedef typename Transformer::result_type result_type;

   result_type operator()(argument_type x) const
   { return transform(x.base(), compose(G(), x.functor())); }

   result_type operator()(argument_type x, G const& g) const
   { return transform(x.base(), compose(g, x.functor())); }
};

template <typename Base, typename F, typename G>
struct TransformMatrix<MatrixTransformProxy<Base, F>, G>
{
   typedef MatrixTransformProxy<Base, F> const& argument_type;
   typedef MatrixTransformProxy<Base, F> const& first_argument_type;
   typedef G second_argument_type;
   typedef Compose<G, F> Composer;
   typedef typename Composer::result_type Func;

   typedef Transform<typename reference_to_arg<Base>::type, Func> Transformer;
   typedef typename Transformer::result_type result_type;

   result_type operator()(argument_type x) const
   { return transform(x.mutable_base(), compose(G(), x.functor())); }

   result_type operator()(argument_type x, G const& g) const
   { return transform(x.mutable_base(), compose(g, x.functor())); }
};

template <typename T, typename F, typename S, typename U>
struct TransformInterface<T, F, Concepts::MatrixExpression<S, U>>
   : TransformMatrix<T, F> {};

template <typename T, typename F, typename S, typename U>
struct TransformRef<T, F, Concepts::MatrixExpression<S, U>>
   : TransformMatrix<T&, F> {};

// if we transform a MatrixTransposeProxy, then shift the transform inside.
// this transformation maybe has a bad effect on compile time.

template <typename Base, typename G>
struct TransformMatrix<MatrixTransposeProxy<Base>, G>
{
   typedef Transform<typename reference_to_arg<Base>::type, G> Transformer;

   BOOST_MPL_ASSERT((is_proxy_reference<typename Transformer::result_type>));

   typedef MatrixTranspose<
      typename reference_to_arg<typename Transformer::result_type>::type
    , Identity<typename make_value<typename G::result_type>::type>
   > Transposer;

   typedef typename Transposer::result_type result_type;
   typedef MatrixTransposeProxy<Base> const& argument_type;
   typedef argument_type first_argument_type;
   typedef G second_argument_type;

   result_type operator()(argument_type m) const
   { return Transposer()(transform(m.mutable_base(), G())); }

   result_type operator()(argument_type m, G const& g) const
   { return Transposer()(transform(m.mutable_base(), g)); }
};

// TODO: lvalue version of TransformMatrix<MatrixTransposeProxy<Base>, G>

// Assign

template <typename LHS, typename Base, typename F>
struct AssignExpression<LHS&, MatrixTransformProxy<Base, F> >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef MatrixTransformProxy<Base, F> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      assign(x, transform(eval_expression(y.base()), y.functor()));
   }
};

// Add

template <typename LHS, typename Base, typename F>
struct AddExpression<LHS&, MatrixTransformProxy<Base, F> >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef MatrixTransformProxy<Base, F> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      add(x, transform(eval_expression(y.base()), y.functor()));
   }
};

// Subtract

template <typename LHS, typename Base, typename F>
struct SubtractExpression<LHS&, MatrixTransformProxy<Base, F> >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef MatrixTransformProxy<Base, F> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      subtract(x, transform(eval_expression(y.base()), y.functor()));
   }
};

// SwapSortOrder

template <typename Base, typename F>
struct SwapSortOrder<MatrixTransformProxy<Base, F> >
{
   typedef typename SwapSortOrder<
      typename reference_to_arg<Base>::type
   >::result_type SwappedBase;

   BOOST_MPL_ASSERT((is_proxy_reference<SwappedBase>));

   typedef typename Transform<
      typename make_const_reference<SwappedBase>::type
    , F
   >::result_type result_type;

  typedef MatrixTransformProxy<Base, F> const& argument_type;
   result_type operator()(argument_type x) const
   { return transform(swap_sort_order(x.base()),x.functor()); }
};

template <typename Base, typename F>
struct SwapSortOrder<MatrixTransformProxy<Base, F>&>
{
   typedef typename SwapSortOrder<
      typename reference_to_arg<Base>::type
   >::result_type SwappedBase;

   BOOST_MPL_ASSERT((is_proxy_reference<SwappedBase>));

   typedef typename Transform<
      typename make_reference<SwappedBase>::type
    , F
   >::result_type result_type;

   typedef MatrixTransformProxy<Base, F>& argument_type;
   result_type operator()(argument_type x) const
   { return transform(swap_sort_order(x.base()),x.functor()); }
};

} // namespace LinearAlgebra

#endif
