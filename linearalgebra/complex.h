// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/complex.h
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
  complex.h

  An 'unwrapped' vector of a complex type, that stores the real and complex parts
  separately.

  The requirements on the types of the real and imag parts is that the
  corresponding Real and Imag functions are identity functions.

  Created 2005-01-13 Ian McCulloch
*/

#if !defined(COMPLEXVECTOR_H_CNKJH78TY983Q5YP3Y89PYQW34P9)
#define COMPLEXVECTOR_H_CNKJH78TY983Q5YP3Y89PYQW34P9

#include "vectoroperationsbase.h"
#include <boost/mpl/assert.hpp>

namespace LinearAlgebra
{

template <typename RealPart, typename ImagPart = RealPart>
class Complex
{
   public:
      typedef RealPart real_type;
      typedef ImagPart imag_type;

      Complex() {}

      Complex(Complex const& c) : Real_(c.Real_), Imag_(c.Imag_) {}

      explicit Complex(real_type const& R) : Real_(R), Imag_() {}

      Complex(real_type const& R, imag_type const& I) : Real_(R), Imag_(I) {}

      template <typename U>
      typename boost::enable_if<is_convertible<U>, real_type>::type 
      explicit Complex(U const& r) : Real_(r) {}

      template <typename U, typename V>
      typename boost::enable_if<boost::mpl::and_<boost::is_convertible<U, real_type>,
						 boost::is_convertible<V, imag_type> > >::type 

      template <typename U>
      typename boost::enable_if<is_convertible<U, real_type> >::type 
      operator=(U const& x)
      {
	 real_type Temp(x);
	 assign(Real_, Temp);
	 zero(Imag_);
	 return *this;
      }

      template <typename U>
      typename boost::enable_if<is_convertible<U, real_type> >::type 
      operator=(NoAliasProxy<U> const& x)
      {
	 assign(Real_, x.value());
	 zero(Imag_);
	 return *this;
      }

      real_type& real() { return Real_; }
      imag_type& imag() { return Imag_; }

      real_type const& real() const { return Real_; }
      imag_type const& imag() const { return Imag_; }

   private:
      real_type Real_;
      imag_type Imag_;
};

// Real

template <typename RealPart, typename ImagPart>
struct Real<Complex<RealPart, ImagPart> >
{
   typedef RealPart const& result_type;
   typedef Complex<RealPart, ImagPart> argument_type;
   result_type operator()(argument_type const& v) const { return v.real(); }
};

// Real&

template <typename RealPart, typename ImagPart>
struct Real<Complex<RealPart, ImagPart>&>
{
   typedef RealPart& result_type;
   typedef Complex<RealPart, ImagPart>& argument_type;
   result_type operator()(argument_type v) const { return v.real(); }
};

// Imag

template <typename RealPart, typename ImagPart>
struct Imag<Complex<RealPart, ImagPart> >
{
   typedef RealPart const& result_type;
   typedef Complex<RealPart, ImagPart> argument_type;
   result_type operator()(argument_type const& v) const { return v.real(); }
};

// Imag&

template <typename RealPart, typename ImagPart>
struct Imag<Complex<RealPart, ImagPart>&>
{
   typedef RealPart& result_type;
   typedef Complex<RealPart, ImagPart>& argument_type;
   result_type operator()(argument_type v) const { return v.real(); }
};

// StreamInsert

template <typename R, typename I>
struct StreamInsert<Complex<R, I> >
{
   typedef std::ostream& result_type;
   typedef std::ostream& first_argument_type;
   typedef Complex<R, I> second_argument_type;
   result_type operator()(std::ostream& out, second_argument_type const& x) const
      { return out << "{real=" << x.real() << ", imag=" << x.imag() << '}'; }
};
   
// assignment

template <typename LHS, typename R, typename I>
struct Assign<LHS&, Complex<R, I> >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef Complex<R, I> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      assign(real(x), y.real());
      assign(imag(x), y.imag());
   }
};

template <typename R, typename I, typename RHS>
struct Assign<Complex<R, I>&, RHS>
{
   typedef void result_type;
   typedef Complex<R, I>& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      assign(x.real(), real(y));
      assign(x.imag(), imag(y));
   }
};

template <typename RL, typename IL, typename RR, typename IR>
struct Assign<Complex<RL, IL>&, Complex<RR, IR> >
{
   typedef void result_type;
   typedef Complex<RL, IL>& first_argument_type;
   typedef Complex<RR, IR> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      assign(x.real(), y.real());
      assign(x.imag(), y.imag());
   }
};

// Add

template <typename LHS, typename R, typename I>
struct Add<LHS&, Complex<R, I> >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef Complex<R, I> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      add(real(x), y.real());
      add(imag(x), y.imag());
   }
};

template <typename R, typename I, typename RHS>
struct Add<Complex<R, I>&, RHS>
{
   typedef void result_type;
   typedef Complex<R, I>& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      add(x.real(), real(y));
      add(x.imag(), imag(y));
   }
};

template <typename RL, typename IL, typename RR, typename IR>
struct Add<Complex<RL, IL>&, Complex<RR, IR> >
{
   typedef void result_type;
   typedef Complex<RL, IL>& first_argument_type;
   typedef Complex<RR, IR> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      add(x.real(), y.real());
      add(x.imag(), y.imag());
   }
};

// Subtract

template <typename LHS, typename R, typename I>
struct Subtract<LHS&, Complex<R, I> >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef Complex<R, I> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      subtract(real(x), y.real());
      subtract(imag(x), y.imag());
   }
};

template <typename R, typename I, typename RHS>
struct Subtract<Complex<R, I>&, RHS>
{
   typedef void result_type;
   typedef Complex<R, I>& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      subtract(x.real(), real(y));
      subtract(x.imag(), imag(y));
   }
};

template <typename RL, typename IL, typename RR, typename IR>
struct Subtract<Complex<RL, IL>&, Complex<RR, IR> >
{
   typedef void result_type;
   typedef Complex<RL, IL>& first_argument_type;
   typedef Complex<RR, IR> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      subtract(x.real(), y.real());
      subtract(x.imag(), y.imag());
   }
};

// Multiply
// This is a bit more subtle, because we need to be careful of aliasing

template <typename LHS, typename R, typename I>
struct Multiply<LHS&, Complex<R, I> >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef Complex<R, I> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      typename make_value<typename Real<LHS&>::result_type>::type Temp = real(x)*real(y) - imag(x)*imag(y);
      imag(x) = noalias(real(x)*imag(y) + imag(x)*real(y));
      real(x) = noalias(Temp);
   }
};

template <typename R, typename I, typename RHS>
struct Multiply<Complex<R, I>&, RHS>
{
   typedef Complex<R, I> LHS;
   typedef void result_type;
   typedef Complex<R, I>& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      typename make_value<typename Real<LHS&>::result_type>::type Temp = real(x)*real(y) - imag(x)*imag(y);
      imag(x) = noalias(real(x)*imag(y) + imag(x)*real(y));
      real(x) = noalias(Temp);
   }
};

template <typename RL, typename IL, typename RR, typename IR>
struct Multiply<Complex<RL, IL>&, Complex<RR, IR> >
{
   typedef void result_type;
   typedef Complex<RL, IL> LHS;
   typedef Complex<RL, IL>& first_argument_type;
   typedef Complex<RR, IR> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      typename make_value<typename Real<LHS&>::result_type>::type Temp = real(x)*real(y) - imag(x)*imag(y);
      imag(x) = real(x)*imag(y) + imag(x)*real(y);
      real(x) = Temp;
   }
};

//
// Addition
//

template <typename LHS, typename R, typename I>
struct Addition<LHS, Complex<R, I> >
{
   typedef typename make_value<typename Real<LHS>::result_type>::type RealType;
   typedef Complex<RealType, RealType> result_type;
   typedef LHS const& first_argument_type;
   typedef Complex<R, I> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      return result_type(real(x)+real(y), imag(x)+imag(y));
   }
};

template <typename R, typename I, typename RHS>
struct Addition<Complex<R, I>, RHS>
{
   typedef typename make_value<typename Real<RHS>::result_type>::type RealType;
   typedef Complex<RealType, RealType> result_type;
   typedef Complex<R, I> const& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      return result_type(real(x)+real(y), imag(x)+imag(y));
   }
};

template <typename R1, typename I1, typename R2, typename I2>
struct Addition<Complex<R1, I1>, Complex<R2, I2> >
{
   typedef typename make_value<typename Addition<R1, R2>::result_type>::type RealType;
   typedef typename make_value<typename Addition<I1, I2>::result_type>::type ImagType;
   typedef Complex<RealType, ImagType> result_type;
   typedef Complex<R1, I1> const& first_argument_type;
   typedef Complex<R2, I2> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      return result_type(real(x)+real(y), imag(x)+imag(y));
   }
};


} // namespace LinearAlgebra

#endif
