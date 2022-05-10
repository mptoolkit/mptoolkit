// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/numerics.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
  numerics.h

  Some functions for reliable numerics.

  Created 2004-02-27 Ian McCulloch

  After getting burned with compliers not defining abs(double) causing implicit conversion to abs(int),
  we define here the euclidean norm norm_2(T), for types T = float (does anyone use float anymore?),
  T = double and T = std::complex<double>, and integer types.
  norm_2 always returns a floating point type.  norm_2_sq(T) gives the square of the 2-norm.

  Also defined are conj(x), which gives sensible results for real types as well as complex,
  and herm(x), which reduces to conj(x) for scalar types.

  pythag(x,y) returns sqrt(x*x + y*y), without extraneous loss of precision.

  TODO: possibly, norm_2, norm_2_sq herm and conj should be defined in the global namespace?  Problem
  is ADL will still not find them, because builtin types don't have an associated namespace.  sigh.

  Integer division:
  template <typename T> class div_result, for types T = int, long and long long,
  has the same purpose the standard library types div_t, ldiv_t and lldiv_t respectively.

  For integer division, the sign of the remainder is a function of the sign of the numerator
  and the sign of the denominator.  Thus there are 16 distinct ways to calculate div(num, denom).
  Clearly, some are more useful than others.  We define 4 here.  The naming convention is
  "divX" where X stands for the sign of the remainder, and is one of:
    'q' (remainder is same sign as quotient)
    'p' (remainder is positive)
    'n' (remainder is same sign as numerator)
    'd' (remainder is same sign as denominator)

  divq(T num, T denom) returns the quotient and remainder of num divided by denom, where the
  quotient is rounded towards negative infinity if num is positive and towards positive infinity
  if num is negative, or equivalently, rounded towards zero if denom is positive
  and rounded away from zero if denom is negative, or equivalently,
  the remainder has the same sign as (num divided by denom).

  divp(T num, T denom) returns the  quotient and remainder of num divided by denom, where the
  quotient is rounded towards zero if num is positive and rounded away from zero
  if num is negative, or equivalently, rounded towards negative infinity if denom is positive
  and rounded towards positive infinity if denom is negative, or equivalently,
  the remainder is non-negative.

  divn(T num, T denom) returns the quotient and remainder of num divided by denom, where the
  quotient is rounded towards zero, or equivalently, the sign of the remainder is the same
  as the sign of num.  This is equivalent to std::div(num, denom).

  divd(T num, T denom) returns the quotient and remainder of num divided by denom, where the
  quotient is rounded towards negative infinity, or equivalently, the remainder has the same sign as denom.

  2005-01-12: Added a FloatIsInt<T> class, that is designed to use with the boost numeric cast
  library; it asserts that the supplied floating type has zero fractional part.

  2006-07-06: Added kahan_sum, for calculating the sum of a list of numbers to accuracy 2*epsilon,
  rather than n*epsilon.
  See http://docs.sun.com/source/806-3568/ncg_goldberg.html#1262
*/

#if !defined(NUMERICS_H_HCIUHU48R4Y89Y9832YR982Y)
#define NUMERICS_H_HCIUHU48R4Y89Y9832YR982Y

#include <complex>
#include <cmath>
#include <stdlib.h>
#include <functional>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_arithmetic.hpp>

namespace numerics
{
inline
double pythag(double x, double y)
{
   double Absx = fabs(x);
   double Absy = fabs(y);
   if (Absx > Absy)
   {
      double Ratio = Absy / Absx;
      return Absx * sqrt(1.0 + Ratio * Ratio);
   }
   else if (Absy == 0.0) return 0.0;
   // else
   double Ratio = Absx / Absy;
   return Absy * sqrt(1.0 + Ratio * Ratio);
}

template <typename FwdIter>
inline
typename std::iterator_traits<FwdIter>::value_type
kahan_sum(FwdIter start, FwdIter end)
{
   typename std::iterator_traits<FwdIter>::value_type s = 0, c = 0;
   while (start != end)
   {
      typename std::iterator_traits<FwdIter>::value_type y = *start - c;
      typename std::iterator_traits<FwdIter>::value_type t = s+y;
      c = (t-s) - y;
      s = t;
      ++start;
   }
   return s;
}

} // namespace numerics

// For complex numbers that have a small imaginary component (smaller than epsilon * real part),
// force it to be zero.
template <typename T>
std::complex<T>
remove_small_imag(std::complex<T> x)
{
   using std::abs;
   if (abs(x.imag()) < 30*std::numeric_limits<T>::epsilon() * abs(x.real()))
      x.imag(0.0);
   return x;
}

//
// norm_2, norm_2_sq
//
// Euclidean norm and square of the euclidean norm
//

inline
double norm_2(int x)
{
   return std::abs(x);
}

inline
double norm_2_sq(int x)
{
   return double(x)*x;
}

inline
double norm_2(long x)
{
   return std::abs(x);
}

inline
double norm_2_sq(long x)
{
   return double(x)*x;
}

#if defined(USE_LONGLONG)
inline
double norm_2(long long x)
{
   return std::abs(double(x));
}
inline
double norm_2_sq(long long x)
{
   return double(x)*x;
}

#endif

inline
float norm_2(float x)
{
   return std::fabs(x);
}

inline
float norm_2_sq(float x)
{
   return x*x;
}

inline
double norm_2(double x)
{
   return std::fabs(x);
}

inline
double norm_2_sq(double x)
{
   return x*x;
}

inline
double norm_2(std::complex<double> x)
{
   return numerics::pythag(x.real(), x.imag());
}

inline
double norm_2_sq(std::complex<double> x)
{
   return norm_2_sq(x.real()) + norm_2_sq(x.imag());
}

//
// norm_inf
//

inline
int norm_inf(int x)
{
   return std::abs(x);
}

inline
long norm_inf(long x)
{
   return std::abs(x);
}

inline
float norm_inf(float x)
{
   return std::abs(x);
}

inline
double norm_inf(double x)
{
   return std::abs(x);
}

inline
double norm_inf(std::complex<double> x)
{
   return numerics::pythag(x.real(), x.imag());
}

//
// conj
//
// Returns the conjugate of a real or complex number
//

inline
float conj(float x)
{
   return x;
}

inline
double conj(double x)
{
   return x;
}

inline
std::complex<double> conj(std::complex<double> x)
{
   return std::conj(x);
}

//
// herm
//
// Returns the hermitian conjugate of a real or complex number
//

template
<typename T>
inline
typename boost::enable_if<boost::is_arithmetic<T>, T>::type
herm(T x)
{
   return x;
}

inline
std::complex<float> herm(std::complex<float> x)
{
   return std::conj(x);
}

inline
std::complex<double> herm(std::complex<double> x)
{
   return std::conj(x);
}

//
// transpose
//
// Returns the transpose of a real or complex number - ie. a null-operation
//

template
<typename T>
inline
typename boost::enable_if<boost::is_arithmetic<T>, T>::type
transpose(T x)
{
   return x;
}

inline
std::complex<double> transpose(std::complex<double> x)
{
   return x;
}



namespace numerics
{

using ::norm_2;
using ::norm_2_sq;
using ::norm_inf;
using ::transpose;
using ::conj;
using ::herm;

//
// a functor to test whether the euclidean norm is smaller than some number v.
//

template <typename T>
struct is_norm_2_smaller_than : public std::unary_function<bool, T>
{
   double value;
   is_norm_2_smaller_than(double v) : value(v) {}

   bool operator()(T const& x) const { return norm_2(x) < value; }
};

//
// complex_ref
//

template <typename T>
class complex_ref
{
   public:
      complex_ref(double& re_, double& im_) : re(re_), im(im_) {}

      complex_ref operator=(std::complex<T> const& x)
      { re = x.real(); im = x.imag(); }

      complex_ref operator=(complex_ref const& x)
      { re = x.re; im = x.im; }

      operator std::complex<T>() const { return std::complex<T>(re, im); }

      T const& real() const { return re; }
      T const& imag() const { return im; }

      T& real() { return re; }
      T& imag() { return im; }

   private:
      T& re;
      T& im;
};

// div

template <typename T>
struct div_result
{
   T quot;
   T rem;

   div_result() {}
   div_result(T const& q, T const& r) : quot(q), rem(r) {}
};

template <>
struct div_result<int>
{
   int quot;
   int rem;

   div_result() {}
   div_result(int q, int r) : quot(q), rem(r) {}
   div_result(std::div_t const& d) : quot(d.quot), rem(d.rem) {}
};

template <>
struct div_result<long>
{
   long quot;
   long rem;

   div_result() {}
   div_result(long q, long r) : quot(q), rem(r) {}
   div_result(ldiv_t const& d) : quot(d.quot), rem(d.rem) {}
};

#if defined(USE_LONGLONG)
template <>
struct div_result<long long>
{
   long long quot;
   long long rem;

   div_result() {}
   div_result(long long q, long long r) : quot(q), rem(r) {}
   div_result(lldiv_t const& d) : quot(d.quot), rem(d.rem) {}
};
#endif

namespace Private
{

template <int pn = 3%(-2), int np = (-3)%2, int nn = (-3)%(-2)>
struct div_helper
{
   template <typename T>
   static div_result<T> divp(T n, T d)
   {
      div_result<T> Res(n/d, n%d);
      if (Res.rem < 0)
      {
         if (d < 0)
         {
            ++Res.quot;
            Res.rem -= d;
         }
         else
         {
            --Res.quot;
            Res.rem += d;
         }
      }
      return Res;
   }

   template <typename T>
   static div_result<T> divq(T n, T d)
   {
      div_result<T> Res(n/d, n%d);
      if (n < 0)
      {
         if ((Res.rem > 0 && d > 0) || (Res.rem < 0 && d <= 0))
         {
            ++Res.quot;
            Res.rem -= d;
         }
      }
      else
      {
         if ((Res.rem > 0 && d < 0) || (Res.rem < 0 && d >= 0))
         {
            --Res.quot;
            Res.rem += d;
         }
      }
      return Res;
   }

   template <typename T>
   static div_result<T> divd(T n, T d)
   {
      div_result<T> Res(n/d, n%d);
      if ((Res.rem < 0 && d >= 0) || (Res.rem > 0 && d < 0))
      {
         --Res.quot;
         Res.rem += d;
      }
      return Res;
   }

   template <typename T>
   static div_result<T> divn(T n, T d)
   {
      return div_result<T>(n/d, n%d);
   }
};

// optimize for the case where the sign of the remainder is the
// same as the sign of the numerator.  This is the case on linux/x86 at least.
// It doesn't actually help much, except for divq() and divn()
template <>
struct div_helper<1,-1,-1>
{
   template <typename T>
   static div_result<T> divp(T n, T d)
   {
      div_result<T> Res(n/d, n%d);
      if (Res.rem < 0)
      {
         if (d < 0)
         {
            ++Res.quot;
            Res.rem -= d;
         }
         else
         {
            --Res.quot;
            Res.rem += d;
         }
      }
      return Res;
   }

   template <typename T>
   static div_result<T> divq(T n, T d)
   {
      div_result<T> Res(n/d, n%d);
      if (n < 0)
      {
         // for this case, we know that Rem <= 0
         if (Res.rem != 0 && d < 0)
         {
            ++Res.quot;
            Res.rem -= d;
         }
      }
      else
      {
         // we know here that Rem >= 0
         if (Res.rem != 0 && d < 0)
         {
            --Res.quot;
            Res.rem += d;
         }
      }
      return Res;
   }

   template <typename T>
   static div_result<T> divd(T n, T d)
   {
      div_result<T> Res(n/d, n%d);
      if ((Res.rem < 0 && d >= 0) || (Res.rem > 0 && d < 0))
      {
         --Res.quot;
         Res.rem += d;
      }
      return Res;
   }

   template <typename T>
   static div_result<T> divn(T n, T d)
   {
      return div_result<T>(n/d, n%d);
   }

};

} // namespace Private

//
// For arbitrary types, we cannot make any assumptions about
// the sign of the remainder, so fall back to the safe versions
//

template <typename T>
inline
div_result<T> divq(T const& n, T const& d)
{
   return Private::div_helper<0,0,0>::divq(n, d);
}

template <typename T>
inline
div_result<T> divp(T const& n, T const& d)
{
   return Private::div_helper<0,0,0>::divp(n, d);
}

template <typename T>
inline
div_result<T> divn(T const& n, T const& d)
{
   return Private::div_helper<0,0,0>::divn(n, d);
}

template <typename T>
inline
div_result<T> divd(T const& n, T const& d)
{
   return Private::div_helper<0,0,0>::divd(n, d);
}

//
// For builtins, we can optimize on the sign
//

inline
div_result<int> divq(int n, int d)
{
   return Private::div_helper<>::divq(n, d);
}

inline
div_result<int> divp(int n, int d)
{
   return Private::div_helper<>::divp(n, d);
}

inline
div_result<int> divn(int n, int d)
{
   return std::div(n,d);
}

inline
div_result<int> divd(int n, int d)
{
   return Private::div_helper<>::divd(n, d);
}

// long

inline
div_result<long> divq(long n, long d)
{
   return Private::div_helper<>::divq(n,d);
}

inline
div_result<long> divp(long n, long d)
{
   return Private::div_helper<>::divp(n, d);
}

inline
div_result<long> divn(long n, long d)
{
   return ldiv(n, d);
}

inline
div_result<long> divd(long n, long d)
{
   return Private::div_helper<>::divd(n, d);
}

// long long

#if defined(USE_LONGLONG)
inline
div_result<long long> divq(long long n, long long d)
{
   return Private::div_helper<>::divq(n,d);
}

inline
div_result<long long> divp(long long n, long long d)
{
   return Private::div_helper<>::divp(n, d);
}

inline
div_result<long long> divn(long long n, long long d)
{
   return lldiv(n, d);
}

inline
div_result<long long> divd(long long n, long long d)
{
   return Private::div_helper<>::divd(n, d);
}
#endif

//
// FloatIsInt - a FloatToInt policy class for boost numeric conversion library
// that enforces that the given floating type has zero fractional part.

#if 0
template <typename T>
struct FloatIsInt
{
   typedef T source_type;
   typedef T argument_type;

   static source_type nearbyint(T x)
   { using std::trunc; PRECONDITION_EQUAL(trunc(x), x); return x; }
};
#endif

} // namespace numerics

#endif
