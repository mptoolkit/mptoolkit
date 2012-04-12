/* -*- C++ -*- $Id$

  vectortransform.h

  Created 2005-01 Ian McCulloch
*/

#if !defined(VECTORTRANSFORM_H_KDSHJCKYH48Y98)
#define VECTORTRANSFORM_H_KDSHJCKYH48Y98

#include "vectoroperationsbase.h"
#include "noalias.h"
#include "stdvector.h"
#include "vectortransformiterator.h"
#include "crtp_vector.h"
#include <boost/mpl/assert.hpp>
#include <boost/mpl/equal_to.hpp>
#include <boost/mpl/comparison.hpp>
#include <boost/mpl/print.hpp>

namespace LinearAlgebra
{

   // experimental, not used
template <typename F, typename Enable = void>
struct FunctorTraits
{
   typedef typename F::result_type const_reference;
   typedef boost::mpl::false_ is_assignable;
};

template <typename F>
struct FunctorTraits<F,
   typename boost::enable_if
   <
      boost::mpl::and_
      <
         boost::is_reference<F>,
         boost::is_const<boost::remove_reference<F> >
      >
   >::type>
{
   typedef typename F::result_type reference;
   typedef typename boost::remove_reference<reference>::type const& const_reference;
   typedef boost::mpl::true_ is_assignable;
};

template <typename BaseProxyReference, typename F>
class VectorTransformProxy : public VectorBase<VectorTransformProxy<BaseProxyReference, F> >
{
   public:
      BOOST_MPL_ASSERT((is_proxy_reference<BaseProxyReference>));

      BOOST_MPL_ASSERT((implies<
			is_mutable_proxy_reference<typename F::argument_type>,
			is_mutable_proxy_reference<BaseProxyReference> >));

      typedef typename F::result_type reference;
      typedef typename boost::mpl::if_<is_proxy_reference<reference>,  
				       make_const_reference<reference>,
				       make_value<reference> >::type::type const_reference;
   //      typedef typename make_const_reference<reference>::type const_reference;
      typedef typename make_value<const_reference>::type value_type;

      typedef is_mutable_proxy_reference<reference> proxy;
      typedef boost::mpl::not_<proxy> const_proxy;
   //      typedef is_const_proxy_reference<reference>   const_proxy;

      typedef BaseProxyReference base_reference;
      typedef typename make_const_reference<base_reference>::type base_const_reference;

       // typedef VectorTransformProxy<base_const_reference, F> const_type;
       // here, because F may be expecting to take its argument type by non-const reference,
       // so we must use base_reference rather than base_const_reference.
       // Instead, we must use
       // typedef VectorTransformProxy const const_type;
       // but this is the default, so we omit it.
       // Alternatively, we could probe the argument_type of F to see if it is safe to pass
       // a const reference.
       // MatrixTransformProxy has exactly the same issue.

      typedef F functor_type;

      // declare the abstract interface; only used if BaseType is an expression
      typedef typename abstract_interface<typename make_value<BaseProxyReference>::type>::type
         abstract_interface;

      using VectorBase<VectorTransformProxy>::operator[];

      VectorTransformProxy(base_reference Base, functor_type const& Func)
	: Base_(Base), Func_(Func) { }

      explicit VectorTransformProxy(base_reference Base)
	: Base_(Base), Func_() { }

      // this conversion ctor handles non-const to const conversions.
      template <typename OtherBase>
      VectorTransformProxy(VectorTransformProxy<OtherBase, F> const& e)
         : Base_(e.base()), Func_(e.functor()) { }

      const_reference operator[](size_type n) const
	 { return Func_(Base_[n]); }

      reference operator[](size_type n)
	 { return Func_(Base_[n]); }

      template <typename U>
      typename boost::enable_if<is_vector<U>, VectorTransformProxy&>::type
      operator=(U const& x)
      {
	 CHECK_EQUAL(this->size(), x.size());
	 // TODO: better temp vector type

	 std::vector<value_type> Temp(x.size());
	 make_vec(Temp) = x;

	 assign(*this, Temp);
	 return *this;
      }

      template <typename U>
      typename boost::enable_if<is_vector<U>, VectorTransformProxy&>::type
      operator=(NoAliasProxy<U> const& x)
      {
	 CHECK_EQUAL(this->size(), x.size());
	 assign(*this, x.value());
	 return *this;
      }

      base_const_reference base() const { return Base_; }
      base_reference base() { return Base_; }

      functor_type const& functor() const { return Func_; }

      base_reference mutable_base() const { return Base_; }

   private:
      base_reference Base_;
      functor_type Func_;
};

namespace Private
{

// compute the interface type of a VectorTransform.  We mostly
// take over the same interface as the base type, but we cannot
// be a StrideVector (except perhaps in some special cases that
// are handled elsewhere by overloads).

template <typename T, typename OldInterface>
struct VectorTransformInterface
{
   typedef T type;
};

template <typename T, typename V, typename U>
struct VectorTransformInterface<T, VECTOR_EXPRESSION(V, U) >
{
   typedef typename T::value_type value_type;
   typedef VECTOR_EXPRESSION(value_type, void) type;
};

template <typename T, typename V, typename U>
struct VectorTransformInterface<T, LOCAL_VECTOR(V, U) >
{
   typedef typename T::value_type value_type;
   typedef LOCAL_VECTOR(value_type, void) type;
};

template <typename T, typename V, typename U>
struct VectorTransformInterface<T, DENSE_VECTOR(V, U) >
{
   typedef typename T::value_type value_type;
   typedef DENSE_VECTOR(value_type, void) type;
};

template <typename T, typename V, typename U>
struct VectorTransformInterface<T, COMPRESSED_VECTOR(V, U) >
{
   typedef typename T::value_type value_type;
   typedef COMPRESSED_VECTOR(value_type, void) type;
};

template <typename T, typename V, typename U>
struct VectorTransformInterface<T, INJECTIVE_VECTOR(V, U) >
{
   typedef typename T::value_type value_type;
   typedef INJECTIVE_VECTOR(value_type, void) type;
};

template <typename T, typename V, typename U>
struct VectorTransformInterface<T, ORDERED_VECTOR(V, U) >
{
   typedef typename T::value_type value_type;
   typedef ORDERED_VECTOR(value_type, void) type;
};

} // namespace Private

// interface

template <typename Base, typename F>
struct interface<VectorTransformProxy<Base, F> >
{
   typedef typename Private::VectorTransformInterface<VectorTransformProxy<Base, F>,
					     typename interface<Base>::type>::type type;
   typedef typename VectorTransformProxy<Base, F>::value_type value_type;
};

// size

template <typename Base, typename F>
struct Size<VectorTransformProxy<Base, F> >
{
   typedef size_type result_type;
   typedef VectorTransformProxy<Base, F> const& argument_type;
   result_type operator()(argument_type x) const
   {
      return size(x.base());
   }
};

// for a VectorTransformProxy of a VECTOR_EXPRESSION, we need to
// supply EvalExpression

template <typename Base, typename F, typename Value>
struct EvalExpression<VectorTransformProxy<Base, F>, 
	    VECTOR_EXPRESSION(Value, 
			     TEMPLATE2(VectorTransformProxy<Base, F>))>
{
   typedef typename EvalExpression<Base>::result_type BaseValueType;
   typedef Transform<BaseValueType, F> Transformer;

   //   typedef typename Transformer::result_type result_type;
   typedef typename make_value<typename Transformer::result_type>::type result_type;
   typedef VectorTransformProxy<Base, F> argument_type;

   result_type operator()(argument_type const& x) const
   {
      return result_type(transform(eval_expression(x.base()), x.functor()));
   }
};

template <typename Base, typename F, typename Value>
struct EvalExpression<VectorTransformProxy<Base, F> const, 
	    VECTOR_EXPRESSION(Value, 
			      TEMPLATE2(VectorTransformProxy<Base, F> const)) >
{
   typedef typename EvalExpression<Base>::result_type BaseValueType;
   typedef Transform<BaseValueType, F> Transformer;

   typedef typename make_value<typename Transformer::result_type>::type result_type;
   typedef VectorTransformProxy<Base, F> argument_type;

   result_type operator()(argument_type const& x) const
   {
      return result_type(transform(eval_expression(x.base()), x.functor()));
   }
};

// iterators

template <typename Base, typename F>
struct Iterate<VectorTransformProxy<Base, F>&>
{
   typedef VectorTransformProxy<Base, F>& argument_type;
   typedef typename iterator<typename basic_type<Base>::type>::type base_iterator;
   typedef VectorTransformIterator<base_iterator, F> result_type;
   result_type operator()(argument_type x) const
   {
      return result_type(iterate(x.base()), x.functor());
   }
};

template <typename Base, typename F>
struct Iterate<VectorTransformProxy<Base, F> >
{
   typedef VectorTransformProxy<Base, F> const& argument_type;
   typedef typename Iterate<typename reference_to_arg<Base>::type>::result_type base_iterator;
   typedef VectorTransformIterator<base_iterator, F> result_type;
   result_type operator()(argument_type x) const
   {
      return result_type(iterate(x.mutable_base()), x.functor());
   }
};

#if 0
// get_element

template <typename T, typename F>
struct GetVectorElement<VectorTransformProxy<T, F> >
{
   typedef typename VectorTransformProxy<T, F>::const_reference result_type;
   typedef VectorTransformProxy<T, F> const& first_argument_type;
   typedef size_type second_argument_type;
   result_type operator()(first_argument_type v, second_argument_type n) const
#endif


template <typename T, typename F>
struct TransformVector
{
   typedef T const& argument_type;
   typedef T const& first_argument_type;
   typedef F second_argument_type;
   typedef VectorTransformProxy<typename make_const_reference<T>::type, F> result_type;


   //   typedef typename boost::mpl::print<T>::type dummy;


   result_type operator()(argument_type E) const 
      { return result_type(E, F()); }

   result_type operator()(first_argument_type E, second_argument_type const& f) const 
      { return result_type(E, f); }
};

template <typename T, typename F>
struct TransformVector<T&, F>
{
   typedef T& argument_type;
   typedef T& first_argument_type;
   typedef F second_argument_type;
   typedef VectorTransformProxy<typename make_reference<T>::type, F> result_type;


   //   typedef typename boost::mpl::print<T>::type dummy;


   result_type operator()(argument_type E) const 
      { return result_type(E, F()); }

   result_type operator()(first_argument_type E, second_argument_type const& f) const 
      { return result_type(E, f); }
};

template <typename Base, typename F, typename G>
struct TransformVector<VectorTransformProxy<Base, F>&, G>
{
   typedef VectorTransformProxy<Base, F>& argument_type;
   typedef VectorTransformProxy<Base, F>& first_argument_type;
   typedef G second_argument_type;
   typedef Compose<G, F> Composer;
   typedef typename Composer::result_type Func;
   
   typedef Transform<typename reference_to_arg<Base>::type, Func> Transformer;
   typedef typename Transformer::result_type result_type;

   //   typedef typename boost::mpl::print<result_type>::type dummy;

   result_type operator()(argument_type x) const 
   { return Transformer()(x.base(), compose(G(), x.functor())); }

   result_type operator()(argument_type x, G const& g) const 
   { return Transformer()(x.base(), compose(g, x.functor())); }
};

template <typename Base, typename F, typename G>
struct TransformVector<VectorTransformProxy<Base, F>, G>
{
   typedef VectorTransformProxy<Base, F> argument_type;
   typedef VectorTransformProxy<Base, F> first_argument_type;
   typedef G second_argument_type;
   typedef Compose<G, F> Composer;
   typedef typename Composer::result_type Func;
   
   typedef Transform<typename basic_type<Base>::type, Func> Transformer;
   typedef typename Transformer::result_type result_type;

   //   typedef typename boost::mpl::print<result_type>::type dummy;

   result_type operator()(argument_type const& x) const 
   { return Transformer()(x.base(), compose(G(), x.functor())); }

   result_type operator()(argument_type const& x, G const& g) const 
   { return Transformer()(x.base(), compose(g, x.functor())); }
};

template <typename T, typename F, typename S, typename U>
struct TransformInterface<T, F, VECTOR_EXPRESSION(S, U)>
   : TransformVector<T, F> {};

// expression assignment

template <typename LHS, typename T, typename U>
struct AssignExpression<LHS&, VectorTransformProxy<T, U> >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef VectorTransformProxy<T, U> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      assign(x, transform(eval_expression(y.base()), y.functor()));
   }
};

template <typename LHS, typename T, typename U>
struct AddExpression<LHS&, VectorTransformProxy<T, U> >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef VectorTransformProxy<T, U> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      add(x, transform(eval_expression(y.base()), y.functor()));
   }
};

template <typename LHS, typename T, typename U>
struct SubtractExpression<LHS&, VectorTransformProxy<T, U> >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef VectorTransformProxy<T, U> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      subtract(x, transform(eval_expression(y.base()), y.functor()));
   }
};

// expression specializations for Negate

template <typename LHS, typename T, typename U>
struct AddExpression<LHS&, VectorTransformProxy<T, Negate<U> > >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef VectorTransformProxy<T, Negate<U> > const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      subtract(x, y.base());
   }
};

template <typename LHS, typename T, typename U>
struct SubtractExpression<LHS&, VectorTransformProxy<T, Negate<U> > >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef VectorTransformProxy<T, Negate<U> > const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      add(x, y.base());
   }
};

} // namespace LinearAlgebra

#endif
