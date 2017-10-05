// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/complexvector.h
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
  complexvector.h

  An 'unwrapped' vector of a complex type, that stores the real and compex parts
  separately.

  The requirements on the types of the real and imag parts is that the
  corresponding Real and Imag functions are identity functions.

  Created 2005-01-13 Ian McCulloch

  This is mostly an experiment to see how easy it is to define a type
  that is simply a 'vector expression' without conforming to any other refined concepts.

  This could easily be generalized to a Complex<anything>, not just a vector.
*/

#if !defined(COMPLEXVECTOR_H_CNKJH78TY983Q5YP3Y89PYQW34P9)
#define COMPLEXVECTOR_H_CNKJH78TY983Q5YP3Y89PYQW34P9

#include "vectoroperationsbase.h"
#include <boost/mpl/assert.hpp>

namespace LinearAlgebra
{

// a metafunction to produce the equivalent complex type
// given the real & imag types.  Wouldn't it be nice if std::complex
// was more generic?
template <typename RealPart, typename ImagPart,
   typename RealInterface = typename interface<RealPart>::type,
   typename ImagInterface = typename interface<ImagPart>::type>
struct ComplexType
{
};

template <>
struct ComplexType<double, double, AnyScalar<double>, AnyScalar<double> >
{
   typedef std::complex<double> type;
};

template <>
struct ComplexType<float, float, AnyScalar<float>, AnyScalar<float> >
{
   typedef std::complex<float> type;
};

template <typename RealPart, typename ImagPart = RealPart>
class ComplexVector
{
   public:
      typedef RealPart real_type;
      typedef ImagPart imag_type;

      typedef typename interface<real_type>::value_type real_value_type;
      typedef typename interface<imag_type>::value_type imag_value_type;

      BOOST_MPL_ASSERT_MSG(is_identity<Real<real_value_type> >::value,
                           ComplexVector_must_have_non_complex_real_part,
                           (RealPart));

      typedef typename make_reference<RealPart>::type real_reference;
      typedef typename make_reference<ImagPart>::type imag_reference;

      typedef typename make_const_reference<RealPart>::type const_real_reference;
      typedef typename make_const_reference<ImagPart>::type const_imag_reference;

      typedef typename ComplexType<real_value_type, imag_value_type>::type value_type;

      ComplexVector(real_reference Real,
                    imag_reference Imag) : Real_(Real), Imag_(Imag) {}

      const_real_reference real() const { return Real_; }
      real_reference real() { return Real_; }

      const_imag_reference imag() const { return Imag_; }
      imag_reference imag() { return Imag_; }

      template <typename U>
      typename boost::enable_if<is_vector<U>, ComplexVector& >::type
      operator=(U const& x)
      {
         // TODO: the vector type should be something better here!
         std::vector<value_type> Temp(x.size());
         assign(Temp, x);
         assign(*this, Temp);
         return *this;
      }

      template <typename U>
      typename boost::enable_if<is_vector<U>, ComplexVector& >::type
      operator=(NoAliasProxy<U> const& x)
      {
         assign(*this, x.value());
         return *this;
      }

      size_type size() const
   { //PRECONDITION_EQUAL(Size<real_type>()(Real_), Size<imag_type>()(Imag_));
        return Size<real_type>()(Real_); }

   private:
      real_type Real_;
      real_type Imag_;
};

// TODO: we need to do something better with make_value for complex types; currently
// it would give a Vector<complex<T> >.
// Perhaps make_value is OK as it is, but we want to do something better for
// the temporary type used in (alias) assignment.

template <typename R, typename I>
struct make_value<ComplexVector<R, I> >
{
   typedef typename ComplexType<R, I>::type type;
};

// generate complex vector from real and imaginary parts

template <typename R, typename I>
typename boost::enable_if<
   boost::mpl::and_<is_vector<R>, is_vector<I> >,
   ComplexVector<typename make_const_reference<R>::type, typename make_const_reference<I>::type>
>::type
make_complex(R const& r, I const& i)
{
   return ComplexVector<typename make_const_reference<R>::type, typename make_const_reference<I>::type>(r,i);
}

template <typename R, typename I>
typename boost::enable_if<
   boost::mpl::and_<is_vector<R>, is_vector<I> >,
   ComplexVector<typename make_reference<R>::type, typename make_reference<I>::type>
>::type
make_complex(R& r, I& i)
{
   return ComplexVector<typename make_reference<R>::type, typename make_reference<I>::type>(r,i);
}

// interface

// TODO: ComplexVector of dense vectors could satisfy the DENSE_VECTOR interface,
// but we would probably still want to specialize assignment for complex left-hand side.

template <typename RealPart, typename ImagPart>
struct interface<ComplexVector<RealPart, ImagPart> >
{
   typedef ComplexVector<RealPart, ImagPart> vector_type;
   typedef typename vector_type::value_type value_type;
   typedef VECTOR_EXPRESSION(value_type, vector_type) type;
};

template <typename RealPart, typename ImagPart>
struct abstract_interface<ComplexVector<RealPart, ImagPart> >
   : vector_abstract_or<typename abstract_interface<RealPart>::type,
                        typename abstract_interface<ImagPart>::type>
{
};

// Real

template <typename RealPart, typename ImagPart, typename V, typename Vi>
struct RealInterface<ComplexVector<RealPart, ImagPart>, VECTOR_EXPRESSION(V, Vi)>
{
   typedef RealPart const& result_type;
   typedef ComplexVector<RealPart, ImagPart> argument_type;
   result_type operator()(argument_type const& v) const { return v.real(); }
};

// Real&

template <typename RealPart, typename ImagPart, typename V, typename Vi>
struct RealInterface<ComplexVector<RealPart, ImagPart>&, VECTOR_EXPRESSION(V, Vi)>
{
   typedef RealPart& result_type;
   typedef ComplexVector<RealPart, ImagPart>& argument_type;
   result_type operator()(argument_type v) const { return v.real(); }
};

// Imag

template <typename RealPart, typename ImagPart, typename V, typename Vi>
struct ImagInterface<ComplexVector<RealPart, ImagPart>, VECTOR_EXPRESSION(V, Vi)>
{
   typedef RealPart const& result_type;
   typedef ComplexVector<RealPart, ImagPart> argument_type;
   result_type operator()(argument_type const& v) const { return v.real(); }
};

// Imag&

template <typename RealPart, typename ImagPart, typename V, typename Vi>
struct ImagInterface<ComplexVector<RealPart, ImagPart>&, VECTOR_EXPRESSION(V, Vi)>
{
   typedef RealPart& result_type;
   typedef ComplexVector<RealPart, ImagPart>& argument_type;
   result_type operator()(argument_type v) const { return v.real(); }
};

// a ComplexVector of T,U has ComplexType<T,U> as its value_type

template <typename R, typename I, typename Rv, typename Ri, typename Iv, typename Ii>
struct ComplexType<R, I, ANY_VECTOR(Rv, Ri), ANY_VECTOR(Iv, Ii)>
{
   typedef typename ComplexType<Rv, Iv>::type value_type;
   typedef typename abstract_interface<ComplexVector<R, I> >::type abstract_type;
   typedef typename make_vector_from_abstract<value_type, abstract_type>::type type;
};

// StreamInsert

template <typename R, typename I, typename V, typename Vi>
struct StreamInsert<ComplexVector<R, I>, VECTOR_EXPRESSION(V, Vi)>
{
   typedef std::ostream& result_type;
   typedef std::ostream& first_argument_type;
   typedef ComplexVector<R, I> second_argument_type;
   result_type operator()(std::ostream& out, second_argument_type const& x) const
      { return out << "{real=" << x.real() << ", imag=" << x.imag() << '}'; }
};

// assignment

template <typename LHS, typename R, typename I>
struct Assign<LHS&, ComplexVector<R, I> >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef ComplexVector<R, I> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      assign(real(x), y.real());
      assign(imag(x), y.imag());
   }
};

template <typename R, typename I, typename RHS>
struct Assign<ComplexVector<R, I>&, RHS>
{
   typedef void result_type;
   typedef ComplexVector<R, I>& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      assign(x.real(), real(y));
      assign(x.imag(), imag(y));
   }
};

template <typename RL, typename IL, typename RR, typename IR>
struct Assign<ComplexVector<RL, IL>&, ComplexVector<RR, IR> >
{
   typedef void result_type;
   typedef ComplexVector<RL, IL>& first_argument_type;
   typedef ComplexVector<RR, IR> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      assign(x.real(), y.real());
      assign(x.imag(), y.imag());
   }
};

// Add

template <typename LHS, typename R, typename I>
struct Add<LHS&, ComplexVector<R, I> >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef ComplexVector<R, I> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      add(real(x), y.real());
      add(imag(x), y.imag());
   }
};

template <typename R, typename I, typename RHS>
struct Add<ComplexVector<R, I>&, RHS>
{
   typedef void result_type;
   typedef ComplexVector<R, I>& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      add(x.real(), real(y));
      add(x.imag(), imag(y));
   }
};

template <typename RL, typename IL, typename RR, typename IR>
struct Add<ComplexVector<RL, IL>&, ComplexVector<RR, IR> >
{
   typedef void result_type;
   typedef ComplexVector<RL, IL>& first_argument_type;
   typedef ComplexVector<RR, IR> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      add(x.real(), y.real());
      add(x.imag(), y.imag());
   }
};

// Subtract

template <typename LHS, typename R, typename I>
struct Subtract<LHS&, ComplexVector<R, I> >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef ComplexVector<R, I> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      subtract(real(x), y.real());
      subtract(imag(x), y.imag());
   }
};

template <typename R, typename I, typename RHS>
struct Subtract<ComplexVector<R, I>&, RHS>
{
   typedef void result_type;
   typedef ComplexVector<R, I>& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      subtract(x.real(), real(y));
      subtract(x.imag(), imag(y));
   }
};

template <typename RL, typename IL, typename RR, typename IR>
struct Subtract<ComplexVector<RL, IL>&, ComplexVector<RR, IR> >
{
   typedef void result_type;
   typedef ComplexVector<RL, IL>& first_argument_type;
   typedef ComplexVector<RR, IR> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      subtract(x.real(), y.real());
      subtract(x.imag(), y.imag());
   }
};

} // namespace LinearAlgebra

#endif
