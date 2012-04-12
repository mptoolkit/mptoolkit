/* -*- C++ -*- $Id$

  fixedvector.h

  A vector type that contains only one element, repeated.

  Created 2005-02-24 Ian McCulloch

  TODO: some more operations could be optimized for this type.  eg, norm.
*/

#if !defined(FIXEDVECTOR_H_HFUIHFU3HYT89YT589Y)
#define FIXEDVECTOR_H_HFUIHFU3HYT89YT589Y

#include "vectoroperations.h"
#include "vectortransform.h"

namespace LinearAlgebra
{

// an iterator that represents a constant sequence
template <typename T>
class ConstantIterator
{
   public:
      typedef T value_type;
      typedef typename make_const_reference<T>::type reference;
      typedef T const* pointer;
      typedef vector_iterator_dense category;

      ConstantIterator() {}

      ConstantIterator(size_type sz, value_type const& x, size_type loc = 0) 
	 : size_(sz), loc_(loc), value_(x) {}

      size_type index() const { return loc_; }

      operator bool() const { return loc_ != size_; }

      pointer operator->() const { DEBUG_CHECK(loc_ != size_); return &value_; }
      reference operator*() const { DEBUG_CHECK(loc_ != size_); return value_; }

      ConstantIterator& operator++() { ++loc_; return *this; }
      ConstantIterator& operator++(int) { return ConstantIterator(size_, value_, loc_++); }

      ConstantIterator& operator+=(size_type n) { loc_ += n; return *this; }

      reference operator[](difference_type n) const { return value_; }

   private:
      size_type size_, loc_;
      value_type value_;
};

template <typename T>
class FixedVector
{
   public:
      typedef is_mutable_proxy_reference<T> proxy;
      typedef is_const_proxy_reference<T> const_proxy;

      typedef T data_type;
      typedef typename boost::add_reference<T const>::type call_type;
      typedef typename make_value<T>::type value_type;
      typedef typename make_reference<T>::type reference;
      typedef typename make_const_reference<T>::type const_reference;

      FixedVector() : size_(0) {}

      FixedVector(size_type sz, call_type val = data_type())
         : size_(sz), value_(val) { }
      
      template <typename U>
      FixedVector(FixedVector<U> const& x)
	 : size_(x.size()), value_(x.value()) {}

      template <typename U>
      FixedVector& operator=(FixedVector<U> const& x)
      {
	 if (proxy::value)
	 {
	    PRECONDITION_EQUAL(size_, x.size());
	 }
	 else
	 {
	    size_ = x.size();
	 }
	 value_ = x.value();
	 return *this;
      }

      size_type size() const { return size_; }

      reference value() { return value_; }
      const_reference value() const { return value_; }

      void resize(size_type newsize) { size_ = newsize; }

   private:
      size_type size_;
      data_type value_;
};

// interface

template <typename T>
struct interface<FixedVector<T> >
{
   typedef typename make_value<T>::type value_type;
   typedef DENSE_VECTOR(value_type, FixedVector<T>) type;
};

// iterators

template <typename T>
struct Iterate<FixedVector<T> >
{
   typedef ConstantIterator<typename make_value<T>::type> result_type;
   typedef FixedVector<T> const& argument_type;
   result_type operator()(argument_type v) const
   {
      return result_type(v.size(), v.value());
   }
};

// resize

template <typename T>
struct Resize<FixedVector<T>&>
{
   typedef void result_type;
   typedef FixedVector<T>& first_argument_type;
   typedef size_type second_argument_type;
   void operator()(FixedVector<T>& v, size_type n) const
   {
      v.resize(n);
   }
};

// GetVectorElement

template <typename T>
struct GetVectorElement<FixedVector<T> >
{
   typedef typename FixedVector<T>::const_reference result_type;
   typedef FixedVector<T> const& first_argument_type;
   typedef size_type second_argument_type;
   result_type operator()(first_argument_type v, second_argument_type i) const
   {
      DEBUG_CHECK(i < v.size());
      return v.value();
   }
};

// Multiply

template <typename T, typename RHS, typename Enable = void>
struct MultiplyFixedVector {};

template <typename T, typename RHS>
struct MultiplyFixedVector<T, RHS, typename boost::enable_if<is_defined<Multiply<T&, RHS> > >::type>
{
   typedef void result_type;
   typedef FixedVector<T>& first_argument_type;
   typedef RHS const& second_argument_type;
   void operator()(FixedVector<T>& v, RHS const& x)
   {
      multiply(v.value(), x);
   }
};

template <typename T, typename RHS, typename RHSi>
struct MultiplyInterface<FixedVector<T>&, RHS, DENSE_VECTOR(T, FixedVector<T>), AnyScalar<RHSi> >
   : MultiplyFixedVector<T, RHS> {};

// Transform

template <typename T, typename F>
struct TransformVector<FixedVector<T>, F>
{
   typedef FixedVector<typename F::result_type> result_type;
   typedef FixedVector<T> const& argument_type;
   typedef FixedVector<T> const& first_argument_type;
   typedef F second_argument_type;

   result_type operator()(argument_type v) const
   { return result_type(v.size(), F()(v.value())); }

   result_type operator()(first_argument_type v, F f) const
   { return result_type(v.size(), f(v.value())); }
};

template <typename T, typename F>
struct TransformVector<FixedVector<T>&, F>
{
   typedef FixedVector<typename F::result_type> result_type;
   typedef FixedVector<T>& argument_type;
   typedef FixedVector<T>& first_argument_type;
   typedef F second_argument_type;

   result_type operator()(argument_type v) const
   { return result_type(v.size(), F()(v.value())); }

   result_type operator()(first_argument_type v, F f) const
   { return result_type(v.size(), f(v.value())); }
};

// BinaryTransform

template <typename S, typename T, typename F, typename Sv, typename Si, typename Tv, typename Ti>
struct BinaryTransform<FixedVector<S>, FixedVector<T>, F, 
                       DENSE_VECTOR(Sv, Si), DENSE_VECTOR(Tv, Ti)>
{
   typedef FixedVector<typename make_value<typename F::result_type>::type> result_type;
   typedef FixedVector<S> const& first_argument_type;
   typedef FixedVector<T> const& second_argument_type;
   typedef F third_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      PRECONDITION_EQUAL(x.size(), y.size());
      return result_type(x.size(), F()(x.value(), y.value()));
   }

   result_type operator()(first_argument_type x, second_argument_type y, F f) const
   {
      PRECONDITION_EQUAL(x.size(), y.size());
      return result_type(x.size(), f(x.value(), y.value()));
   }
};

} // namespace LinearAlgebra

#endif
