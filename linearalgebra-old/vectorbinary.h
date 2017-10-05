// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/vectorbinary.h
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
  vectorbinary.h

  A proxy class & iterators for binary expressions of vectors

  Created 2005-01-10 Ian McCulloch
*/

#if !defined(VECTORBINARY_H_KDSHJCKYH48Y98)
#define VECTORBINARY_H_KDSHJCKYH48Y98

namespace LinearAlgebra
{

template <typename T1, typename T2, typename F>
class VectorBinaryTransformProxy
{
   public:
      // mark this type as a const proxy reference
      typedef boost::mpl::true_ const_proxy;

      // the abstract interface type -
      // shortcut for specializing LinearAlgebra::abstract_interface
      typedef vector_abstract_or<typename abstract_interface<T1>::type,
                                 typename abstract_interface<T2>::type> abstract_interface;

      typedef typename interface<T1>::value_type value1_type;
      typedef typename interface<T2>::value_type value2_type;

      typedef F functor_type;
      typedef typename functor_type::result_type reference;
      typedef typename make_const_reference<reference>::type const_reference;
      typedef typename make_value<reference>::type value_type;
      typedef typename make_const_reference<T1>::type reference1;
      typedef typename make_const_reference<T2>::type reference2;

      VectorBinaryTransformProxy(reference1 x, reference2 y, functor_type f)
         : x_(x), y_(y), f_(f) {} //{ CHECK_EQUAL(Size<T1>()(x_), Size<T2>()(y_)); }

      size_type size() const { using LinearAlgebra::size; return size(x_); }
   // return Size<T1>()(x_); }

      const_reference operator[](size_type n) const { return f_(x_[n], y_[n]); }

      reference1 vector1() const { return x_; }
      reference2 vector2() const { return y_; }

      functor_type const& functor() const { return f_; }

   private:
      reference1 x_;
      reference2 y_;
      functor_type f_;
};

template <typename T1, typename T2, typename Func,
          typename I1 = typename interface<T1>::type,
          typename I2 = typename interface<T2>::type,
          typename Value = typename make_value<typename Func::result_type>::type>
struct VectorBinaryInterface;

template <typename T1, typename T2, typename Func,
          typename I1v, typename I1i,
          typename I2v, typename I2i,
          typename Value>
struct VectorBinaryInterface<T1, T2, Func,
                             DENSE_VECTOR(I1v, I1i),
                             DENSE_VECTOR(I2v, I2i), Value>
{
   typedef DENSE_VECTOR(Value, void) type;
};

template <typename T1, typename T2, typename F>
struct interface<VectorBinaryTransformProxy<T1, T2, F> >
   : public VectorBinaryInterface<T1, T2, F>
{
};

// iterators

template <typename T1, typename T2, typename F>
struct Iterate<VectorBinaryTransformProxy<T1, T2, F> >
{
   typedef typename const_iterator<T1>::type iter1_type;
   typedef typename const_iterator<T2>::type iter2_type;
   typedef VectorBinaryTransformProxy<T1, T2, F> const& argument_type;
   typedef typename VectorBinaryTransformProxy<T1, T2, F>::functor_type func_type;
   typedef VectorBinaryIterator<iter1_type, iter2_type, func_type> result_type;
   result_type operator()(argument_type x) const
   {
      return result_type(iterate(x.vector1()), iterate(x.vector2()), x.functor());
   }
};

template <typename S, typename T, typename F, typename Sv, typename Si, typename Tv, typename Ti>
struct BinaryTransformVector<S, T, F, DENSE_VECTOR(Sv, Si), DENSE_VECTOR(Tv, Ti) >
{
   typedef typename is_semiregular<F>::type semiregular;
   typedef typename is_regular<F>::type regular;
   typedef VectorBinaryTransformProxy<S, T, F> result_type;
   typedef S first_argument_type;
   typedef T second_argument_type;
   typedef F const& third_argument_type;
   result_type operator()(S const& x, T const& y) const { return result_type(x, y); }

   result_type operator()(S const& x, T const& y, F const& f) const
   { return result_type(x, y, f); }
};

template <typename S, typename T, typename F, typename Sv, typename Si, typename Tv, typename Ti>
struct BinaryTransformVectorSemiregular<S, T, F, ORDERED_VECTOR(Sv, Si), ORDERED_VECTOR(Tv, Ti) >
{
   typedef typename is_semiregular<F>::type semiregular;
   typedef typename is_regular<F>::type regular;
   typedef VectorBinaryTransformProxy<S, T, F> result_type;
   typedef S first_argument_type;
   typedef T second_argument_type;
   typedef F const& third_argument_type;
   result_type operator()(S const& x, T const& y) const { return result_type(x, y); }

   result_type operator()(S const& x, T const& y, F const& f) const
   { return result_type(x, y, f); }
};

} // namespace LinearAlgebra

#endif
