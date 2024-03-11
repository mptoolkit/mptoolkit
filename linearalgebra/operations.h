// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/operations.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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
  functions for binary and unary operators & functions.

  Created 2004-04-11 Ian McCuloch

  Now uses a partial specialization based dispatch mechanism.
*/

#if !defined(OPERATIONS_H_CHUIRFHUIYT)
#define OPERATIONS_H_CHUIRFHUIYT

#include "interface.h"
#include "value_with_zero.h"
#include <boost/utility/enable_if.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/mpl/equal_to.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/or.hpp>
#include <boost/mpl/not.hpp>
#include <boost/mpl/and.hpp>
#include <boost/type_traits.hpp>
#include <limits>
#include <cstdlib>
#include <cmath>

// debugging/precondition macro, extension of common/trace.h

#if defined(NDEBUG) || defined(TRACER_DISABLE_ALL_CHECKS)
#define DEBUG_CHECK_CLOSE(x,y)                  \
if (true) ;                                     \
else INVOKE_DUMMY("")
#else
#define DEBUG_CHECK_CLOSE(x, y)                                                 \
if (equal((x), (y))) ;                                                          \
else INVOKE_ASSERT("CHECK_CLOSE", "equal(" #x ", " #y ") is false")(x)(y)
#endif
#if defined(TRACER_DISABLE_ALL_CHECKS)
#define CHECK_CLOSE(x,y)                        \
if (true) ;                                     \
else INVOKE_DUMMY("")
#else
#define CHECK_CLOSE(x, y)                                                       \
if (equal((x), (y))) ;                                                          \
else INVOKE_ASSERT("CHECK_CLOSE", "equal(" #x ", " #y ") is false")(x)(y)
#endif

namespace LinearAlgebra
{

// operator characteristics

// is_builtin.  Used to mark an operator as builtin, so
// don't define it again (eg, for builtin arithmetic operators).  To use, add a nested
// typedef boost::mpl::true_ builtin; to the functor class, or specialize is_builtin.

template <typename T, typename Enable = void>
struct is_builtin : public boost::mpl::false_ { };

template <typename T>
struct is_builtin<T, typename boost::enable_if<typename T::builtin>::type>
   : public boost::mpl::true_ { };

// is_identity.  This is used to eliminate abstraction overhead.
// To mark a functor as an identity functor, add
// typedef boost::mpl::true_ identity; to the functor class.

template <typename T, typename Enable = void>
struct is_identity : public boost::mpl::false_
{
};

template <typename T>
struct is_identity<T, typename boost::enable_if<typename T::identity>::type>
   : public boost::mpl::true_
{
};

// is_stateless - experimental, not used yet

template <typename T, typename Enable = void>
struct is_stateless : public boost::mpl::false_
{
};

template <typename T>
struct is_stateless<T, typename boost::enable_if<typename T::stateless>::type>
   : public boost::mpl::true_
{
};

template <typename T>
struct is_stateless<
   T
 , typename boost::enable_if<
      boost::mpl::or_<
         is_identity<T>
       , is_builtin<T>
      >
   >::type
>
   : boost::mpl::true_ {};

// is_linear.  This marks that a function is linear with respect to addition,
// ie. f(a+b) = f(a) + f(b).  This is important for non-injective
// containers.
// To mark a functor as linear, add
// typedef boost::mpl::true_ linear; to the functor class.

template <typename T, typename Enable = void>
struct is_linear : public boost::mpl::false_
{
};

template <typename T>
struct is_linear<T, typename boost::enable_if<typename T::linear>::type>
   : public boost::mpl::true_
{
};

// is_involutary.  This marks a function satisfies f^2 = identity
// To mark a functor as involutary, add
// typedef boost::mpl::true_ involutary; to the functor class.
// The composition of an involutary function with itself is the identity function.
// The inverse of an involutary function is itself.

template <typename T, typename Enable = void>
struct is_involutary : public boost::mpl::false_
{
};

template <typename T>
struct is_involutary<T, typename boost::enable_if<typename T::involutary>::type>
   : public boost::mpl::true_
{
};

// is_idempotent.  This marks a function satisfies f^2 = f.
// To mark a functor as idempotent, add
// typedef boost::mpl::true_ idempotent; to the functor class.
// The composition of an idempotent function with itself is itself.

template <typename T, typename Enable = void>
struct is_idempotent : public boost::mpl::false_
{
};

template <typename T>
struct is_idempotent<T, typename boost::enable_if<typename T::idempotent>::type>
   : public boost::mpl::true_
{
};

//
// traits for binary functions
//

// is_commutative.
// To mark a binary functor as commutative, add
// typedef boost::mpl::true_ commutative; to the functor class.
// This requires that first_argument_type and second_argument_type are identical.

template <typename T, typename Enable = void>
struct is_commutative : public boost::mpl::false_
{
};

template <typename T>
struct is_commutative<T, typename boost::enable_if<typename T::commutative>::type>
   : public boost::mpl::true_
{
};

// is_bilinear.  This marks that a binary function is bilinear with respect to addition,
// ie. f(a*b,c*d) = a*c*f(b,d).
// To mark a functor as bilinear, add
// typedef boost::mpl::true_ bilinear; to the functor class.

template <typename T, typename Enable = void>
struct is_bilinear : public boost::mpl::false_
{
};

template <typename T>
struct is_bilinear<T, typename boost::enable_if<typename T::bilinear>::type>
   : public boost::mpl::true_
{
};

// is_regular.  This marks that a binary function is 'regular', or
// f(x,0) = f(0,x) = 0.
// To mark a functor as regular, add
// typedef boost::mpl::true_ regular; to the functor class.
// is_bilinear implies is_regular.

template <typename T, typename Enable = void>
struct is_regular : public is_bilinear<T> {};

template <typename T>
struct is_regular<T, typename boost::enable_if<typename T::regular>::type>
   : public boost::mpl::true_ {};

// is_semiregular.  This marks that a binary function is 'semiregular', or
// f(0,0) = 0.
// To mark a functor as semiregular, add
// typedef boost::mpl::true_ semiregular; to the functor class.
// is_regular implies is_semiregular.

template <typename T, typename Enable = void>
struct is_semiregular : public is_regular<T> {};

template <typename T>
struct is_semiregular<T, typename boost::enable_if<typename T::semiregular>::type>
   : public boost::mpl::true_ {};

// is_defined.  true_ if the given function contains
// a result_type, false_ otherwise.

template <typename T, typename Enable = void>
struct is_defined : public boost::mpl::false_ { };

template <typename T>
struct is_defined<T, typename boost::enable_if<exists<typename T::result_type> >::type>
   : public boost::mpl::true_ { };

// is_noalias.  This marks a *container* has assumed to never be aliased.
// This is true for some trivial containers (eg, Range & Slice),
// and also for NoAliasProxy.

template <typename T>
struct is_noalias;

//
// default_tolerance, set_default_tolerance
//
// default_tolerance() returns the tolerance to use for floating point comparisons
//   if the tol parameter is not specified.  This is initially set to
//   std::numeric_limits<double>::epsilon() * 64.
//
// set_default_tolerance(double x) sets the default tolerance to be x.
//

double default_tolerance();

void set_default_tolerance(double x);

// zero

template <typename T, typename tInterface = typename interface<T>::type>
struct ZeroInterface {};

template <typename T, typename Enable>  // primary is in interface.h
struct Zero : ZeroInterface<T> {};

template <typename T>
inline
typename Zero<T>::result_type
zero()
{
   return Zero<T>()();
}

template <typename T, typename Enable>  // primary is in interface.h
struct has_zero : public boost::mpl::false_ {};

template <typename T>
struct has_zero<T, typename boost::enable_if<exists<typename Zero<T>::result_type> >::type>
: boost::mpl::true_ {};

// zero_or_die
// The ordinary zero() function only exists if the type
// has a defined zero value.  However, there are some circumstances
// where we need a zero value, but we don't know until runtime
// whether it will actually be used or not.  Zero_or_die() is
// a version of zero() that always exists, but will abort at runtime
// if there is no zero value available.
// Example: In an inner product, we need to be able to generate a zero
// value, but *only* for the special case where the vectors are empty or do not overlap.
// But it is well defined to take the inner product of vectors
// even if no zero is available, as long as the user makes sure the vectors
// are never empty.

template <typename T, typename Enable = void>
struct ZeroOrDie
{
   typedef T result_type;
   result_type const& operator()() const
   { PANIC("No zero value exists for type")(typeid(T).name());
     return *static_cast<result_type*>(0); }
};

template <typename T>
struct ZeroOrDie<T, typename boost::enable_if<
   exists<typename Zero<T>::result_type> >::type>
   : Zero<T>
{
};

template <typename T>
inline
typename ZeroOrDie<T>::result_type
zero_or_die()
{
   return ZeroOrDie<T>()();
}

// StaticZero
// A statically declared zero value

template <typename T, typename Enable = void>
struct StaticZero { };

template <typename T>
struct StaticZero<T, typename boost::enable_if<exists<typename Zero<T>::result_type> >::type>
{
   static T value;
};

template <typename T>
T
StaticZero<T, typename boost::enable_if<exists<typename Zero<T>::result_type> >::type>
::value = Zero<T>()();

// static_zero_or_die
// returns a const-reference to a statically defined zero value, or
// aborts if no such value exists

template <typename T, typename Enable = void>
struct StaticZeroOrDie
{
   typedef T const& result_type;
   result_type operator()() const
   {
      PANIC("No zero value exists for type")(typeid(T).name());
      abort();
      //return *static_cast<T const*>(0);
   }
};

template <typename T>
struct StaticZeroOrDie<T, typename boost::enable_if<
   exists<typename Zero<T>::result_type> >::type>
{
   typedef T const& result_type;
   result_type operator()() const
   { return StaticZero<T>::value; }
};

template <typename T>
inline
typename StaticZeroOrDie<T>::result_type
static_zero_or_die()
{
   return StaticZeroOrDie<T>()();
}

// zero_or_default
//
// for sparse containers of elements that have no zero, lvalue insertion is
// problematic.  But since it is so useful, we provide a default to use
// if there is no zero value.  The only possible operation on such
// a default is to immediately overwrite it with something else - it could
// be a NaN, or some kind of singular value.

template <typename T, typename tInterface = typename interface<T>::type>
struct DefaultInterface {};

template <typename T, typename Enable = void>
struct Default
{
   typedef T result_type;
   T operator()() const { return T(); }
};

template <typename T>
inline
typename Default<T>::result_type
default_value()
{
   return Default<T>()();
}

template <typename T, typename Enable = void>
struct ZeroOrDefault : Default<T> {};

template <typename T>
struct ZeroOrDefault<T, typename boost::enable_if<
   exists<typename Zero<T>::result_type> >::type>
   : Zero<T>
{
};

template <typename T>
inline
typename ZeroOrDefault<T>::result_type
zero_or_default()
{
   return ZeroOrDefault<T>()();
}

// unary functions

template <typename T>
struct Identity
{
   typedef boost::mpl::true_ identity;

   //typedef typename boost::mpl::print<T>::type dummy;

   typedef T const& result_type;
   typedef T const& argument_type;
   T const& operator()(T const& x) const { return x; }
};

template <typename T>
struct Identity<T&>
{
   typedef boost::mpl::true_ identity;
   typedef T& result_type;
   typedef T& argument_type;
   T& operator()(T& x) const { return x; }
};

struct IdentityF
{
   template <typename T>
   struct result { typedef typename Identity<T>::result_type type; };

   template <typename T>
   typename Identity<T>::result_type operator()(T const& x) const
   {
      return Identity<T>()(x);
   }
};

template <typename T>
inline
typename Identity<T>::result_type
identity(T const& x)
{
   return Identity<T>()(x);
}

template <typename T>
inline
typename Identity<T&>::result_type
identity(T& x)
{
   return Identity<T&>()(x);
}

template <typename T>
struct IdentityRef : Identity<T&> {};

struct IdentityRefF
{
   template <typename T>
   struct result { typedef typename IdentityRef<T>::result_type type; };

   template <typename T>
   typename IdentityRef<T>::result_type operator()(T& x) const
   {
      return Identity<T>()(x);
   }
};

// Inverse

template <typename T, typename TInterface = typename interface<T>::type,
          typename Enable = void>
struct Inverse { };

struct InverseF
{
   template <typename T>
   struct result { typedef typename Inverse<T>::result_type type; };

   template <typename T>
   typename Inverse<T>::result_type operator()(T const& x) const
   {
      return Inverse<T>()(x);
   }
};

template <typename T>
inline
typename T::result_type
inverse(T const& x)
{
   return Inverse<T>(x);
}

// inverse of involutary functions is itself

template <typename T, typename Ti>
struct Inverse<T, Ti, typename boost::enable_if<is_involutary<T> >::type>
{
   typedef T argument_type;
   typedef T result_type;
   result_type operator()(argument_type const& x) const { return x; }
};

// Negate
// In most/all cases, this should be involutary.

template <typename T, typename TInterface = typename interface<T>::type, typename Enable = void>
struct NegateInterface {};

template <typename T,  typename Enable = void>
struct Negate : NegateInterface<T> {};

struct NegateF
{
   template <typename T>
   struct result { typedef typename Negate<T>::result_type type; };

   template <typename T>
   typename Negate<T>::result_type operator()(T const& x) const
   {
      return Negate<T>()(x);
   }
};

template <typename T>
inline
typename boost::disable_if<is_builtin<Negate<T> >,
   typename Negate<T>::result_type>::type
operator-(T const& x)
{
   return Negate<T>()(x);
}

// Trace

template <typename T, typename TInterface = typename interface<T>::type,
          typename Enable = void>
struct TraceInterface
{
};

template <typename T,  typename Enable = void>
struct Trace : TraceInterface<T> {};

struct TraceF
{
   template <typename T>
   struct result { typedef typename Trace<T>::result_type type; };

   template <typename T>
   typename Trace<T>::result_type operator()(T const& x) const
   {
      return Trace<T>()(x);
   }
};

template <typename T>
inline
typename Trace<T>::result_type
trace(T const& x)
{
   return Trace<T>()(x);
}

// Norm2Sq

template <typename T, typename TInterface = typename interface<T>::type,
          typename Enable = void>
struct Norm2Sq
{
};

struct Norm2SqF
{
   template <typename T>
   struct result { typedef typename Norm2Sq<T>::result_type type; };

   template <typename T>
   typename Norm2Sq<T>::result_type operator()(T const& x) const
   {
      return Norm2Sq<T>()(x);
   }
};

template <typename T>
inline
typename Norm2Sq<T>::result_type
norm_2_sq(T const& x)
{
   return Norm2Sq<T>()(x);
}

// Sum

template <typename T, typename TInterface = typename interface<T>::type,
          typename Enable = void>
struct Sum
{
};

struct SumF
{
   template <typename T>
   struct result { typedef typename Sum<T>::result_type type; };

   template <typename T>
   typename Sum<T>::result_type operator()(T const& x) const
   {
      return Sum<T>()(x);
   }
};

template <typename T>
inline
typename Sum<T>::result_type
sum(T const& x)
{
   return Sum<T>()(x);
}

// Norm1 - the 1-norm

template <typename T, typename TInterface = typename interface<T>::type,
          typename Enable = void>
struct Norm1
{
};

struct Norm1F
{
   template <typename T>
   struct result { typedef typename Norm1<T>::result_type type; };

   template <typename T>
   typename Norm1<T>::result_type operator()(T const& x) const
   {
      return Norm1<T>()(x);
   }
};

template <typename T>
inline
typename Norm1<T>::result_type
norm_1(T const& x)
{
   return Norm1<T>()(x);
}

// Norm2 - the 2-norm

// Norm2 has a default implementation, as sqrt(norm_2_sq(x)).
template <typename T, typename TInterface = typename interface<T>::type,
          typename Enable = void>
struct Norm2
{
   typedef typename Norm2Sq<T, TInterface>::result_type result_type;
   typedef T const& argument_type;
   result_type operator()(T x) const { return sqrt(norm_2_sq(x)); }
};

struct Norm2F
{
   template <typename T>
   struct result { typedef typename Norm2<T>::result_type type; };

   template <typename T>
   typename Norm2<T>::result_type operator()(T const& x) const
   {
      return Norm2<T>()(x);
   }
};

template <typename T>
inline
typename Norm2<T>::result_type
norm_2(T const& x)
{
   return Norm2<T>()(x);
}

// NormInf - the infinity-norm

template <typename T, typename TInterface = typename interface<T>::type,
          typename Enable = void>
struct NormInf
{
};

struct NormInfF
{
   template <typename T>
   struct result { typedef typename NormInf<T>::result_type type; };

   template <typename T>
   typename NormInf<T>::result_type operator()(T const& x) const
   {
      return NormInf<T>()(x);
   }
};

template <typename T>
inline
typename NormInf<T>::result_type
norm_inf(T const& x)
{
   return NormInf<T>()(x);
}

// NormFrob - Fronebius norm.
// For vectors and scalars this is equivalent to the 2-norm.
// NormFrob has a default implementation as the square root
// of NormFromSq

template <typename T, typename TInterface = typename interface<T>::type>
struct NormFrobSq {};

struct NormFrobSqF
{
   template <typename T>
   struct result { typedef typename NormFrobSq<T>::result_type type; };

   template <typename T>
   typename NormFrobSq<T>::result_type operator()(T const& x) const
   {
      return NormFrobSq<T>()(x);
   }
};

template <typename T>
inline
typename NormFrobSq<T>::result_type
norm_frob_sq(T const& x)
{
   return NormFrobSq<T>()(x);
}

template <typename T, typename TInterface = typename interface<T>::type,
          typename Enable = void>
struct NormFrobDefault {};

template <typename T, typename TInterface = typename interface<T>::type>
struct NormFrob : NormFrobDefault<T, TInterface> {};

template <typename T, typename TInterface>
struct NormFrobDefault<T, TInterface, typename boost::enable_if<
   is_defined<NormFrobSq<T, TInterface> > >::type>
{
   typedef typename NormFrobSq<T, TInterface>::result_type result_type;
   typedef T const& argument_type;
   result_type operator()(T const& x) const { using std::sqrt; return sqrt(norm_frob_sq(x)); }
};

struct NormFrobF
{
   template <typename T>
   struct result { typedef typename NormFrob<T>::result_type type; };

   template <typename T>
   typename NormFrob<T>::result_type operator()(T const& x) const
   {
      return NormFrob<T>()(x);
   }
};

template <typename T>
inline
typename NormFrob<T>::result_type
norm_frob(T const& x)
{
   return NormFrob<T>()(x);
}

// norm_err
// This is the 'error norm', used to detect cancellation in sparse operations.  It defaults
// to the Frobenius norm.

template <typename T>
decltype(norm_frob_sq(std::declval<T>()))
norm_err_sq(T const& x)
{
   return norm_frob_sq(x);
}

template <typename T>
decltype(norm_err_sq(std::declval<T>()))
norm_err(T const& x)
{
   using std::sqrt;
   return sqrt(norm_err_sq(x));
}

// Transpose

template <typename T, typename TInterface = typename interface<T>::type>
struct TransposeInterface {};

template <typename T>
struct Transpose : TransposeInterface<T> {};

struct TransposeF
{
   template <typename T>
   struct result { typedef typename Transpose<T>::result_type type; };

   template <typename T>
   typename Transpose<T>::result_type operator()(T const& x) const
   {
      return Transpose<T>()(x);
   }
};

template <typename T>
inline
typename Transpose<T>::result_type
transpose(T const& x)
{
   return Transpose<T>()(x);
}

template <typename T>
inline
typename Transpose<T>::result_type
trans(T const& x)
{
   return Transpose<T>()(x);
}

// Conj

template <typename T, typename TInterface = typename interface<T>::type,
          typename Enable = void>
struct ConjInterface {};

template <typename T, typename Enable = void>
struct Conj : ConjInterface<T> {};

struct ConjF
{
   template <typename T>
   struct result { typedef typename Conj<T>::result_type type; };

   template <typename T>
   typename Conj<T>::result_type operator()(T const& x) const
   {
      return Conj<T>()(x);
   }
};

template <typename T>
inline
typename Conj<T>::result_type
conj(T const& x)
{
   return Conj<T>()(x);
}

// The inverse of conj is itself

template <typename T>
struct Inverse<Conj<T> > : public Identity<Conj<T> >
{
};

// Real

template <typename T, typename TInterface = typename interface<T>::type>
struct RealInterface {};

template <typename T>
struct Real : RealInterface<T> {};

template <typename T>
inline
typename Real<typename proxy_lvalue<T>::type>::result_type
real(T const& x)
{
   return Real<typename proxy_lvalue<T>::type>()(proxy_lvalue_cast(x));
}

template <typename T>
inline
typename Real<T&>::result_type
real(T& x)
{
   return Real<T&>()(x);
}

// Imag

template <typename T, typename TInterface = typename interface<T>::type,
          typename Enable = void>
struct ImagInterface {};

template <typename T, typename Enable = void>
struct Imag : ImagInterface<T> {};

template <typename T>
inline
typename Imag<typename proxy_lvalue<T>::type>::result_type
imag(T const& x)
{
   return Imag<typename proxy_lvalue<T>::type>()(proxy_lvalue_cast(x));
}

template <typename T>
inline
typename Imag<T&>::result_type
imag(T& x)
{
   return Imag<T&>()(x);
}

// Herm

template <typename T, typename TInterface = typename interface<T>::type,
          typename Enable = void>
struct HermInterface {};

template <typename T, typename Enable = void>
struct Herm : HermInterface<T> {};

template <typename T>
inline
typename Herm<T>::result_type
herm(T const& x)
{
   return Herm<T>()(x);
}

// The inverse of herm is itself

template <typename T>
struct Inverse<Herm<T> > : public Identity<Herm<T> >
{
};

// Abs

template <typename T, typename TInterface = typename interface<T>::type, typename Enable = void>
struct AbsInterface {};

template <typename T, typename Enable = void>
struct Abs : AbsInterface<T> {};

template <typename T>
inline
typename Abs<T>::result_type
abs(T const& x)
{
   return Abs<T>()(x);
}

// Arg

template <typename T, typename TInterface = typename interface<T>::type, typename Enable = void>
struct ArgInterface {};

template <typename T, typename Enable = void>
struct Arg : ArgInterface<T> {};

template <typename T>
inline
typename Arg<T>::result_type
arg(T const& x)
{
   return Arg<T>()(x);
}

// Sin

template <typename T, typename TInterface = typename interface<T>::type, typename Enable = void>
struct SinInterface {};

template <typename T, typename Enable = void>
struct Sin : SinInterface<T> {};

template <typename T>
inline
typename Sin<T>::result_type
sin(T const& x)
{
   return Sin<T>()(x);
}

// Cos

template <typename T, typename TInterface = typename interface<T>::type, typename Enable = void>
struct CosInterface {};

template <typename T, typename Enable = void>
struct Cos : CosInterface<T> {};

template <typename T>
inline
typename Cos<T>::result_type
cos(T const& x)
{
   return Cos<T>()(x);
}

// Exp

template <typename T, typename TInterface = typename interface<T>::type, typename Enable = void>
struct ExpInterface {};

template <typename T, typename Enable = void>
struct Exp : ExpInterface<T> {};

template <typename T>
inline
typename Exp<T>::result_type
exp(T const& x)
{
   return Exp<T>()(x);
}

// Sqrt

template <typename T, typename TInterface = typename interface<T>::type, typename Enable = void>
struct SqrtInterface {};

template <typename T, typename Enable = void>
struct Sqrt : SqrtInterface<T> {};

template <typename T>
inline
typename Sqrt<T>::result_type
sqrt(T const& x)
{
   return Sqrt<T>()(x);
}

// UnaryComposer

template <typename T1, typename T2, typename Enable = void>
struct UnaryComposer
{
   UnaryComposer(T1 const& f1_ = T1(), T2 const& f2_ = T2()) : first(f1_), second(f2_) {}

   typedef typename T2::argument_type argument_type;
   typedef typename T2::result_type T2Result;

   typedef typename std::result_of<T1(T2Result)>::type result_type;

   result_type operator()(typename T2::argument_type const& x) const
   { return first(second(x)); }

   template <typename T>
   auto operator()(T& x) const -> decltype(first(second(x)))
   { return first(second(x)); }

   T1 first;
   T2 second;
};

// this is a 'new' version.  It won't work because it doesn't have a result_type.
// We need to remove all references to result_type and use std::result_of instead,
// but we also need a method to remove proxies.
template <typename T1, typename T2, typename Enable = void>
struct UnaryComposerX
{
   UnaryComposerX(T1 const& f1_ = T1(), T2 const& f2_ = T2()) : first(f1_), second(f2_) {}

   template <typename T>
   auto operator()(T const& x) const -> decltype(first(second(x))) const
   { return first(second(x)); }

   template <typename T>
   auto operator()(T& x) const -> decltype(first(second(x))) const
   { return first(second(x)); }

   T1 first;
   T2 second;
};

// Swap

template <typename T, typename U,
          typename TInterface = typename interface<T>::type,
          typename UInterface = typename interface<U>::type>
struct SwapInterface {};

template <typename T, typename U, typename Enable = void>
struct Swap : SwapInterface<T, U> {};

template <typename T, typename U>
typename Swap<T&, U&>::result_type
swap(T& t, U& u)
{
   return Swap<T&, U&>()(t,u);
}

template <typename T, typename U>
typename boost::enable_if<is_mutable_proxy<T>, Swap<T&, U&> >::type::result_type
swap(T const& t, U& u)
{
   return Swap<T&, U&>()(const_cast<T&>(t),u);
}

template <typename T, typename U>
typename boost::enable_if<is_mutable_proxy<U>, Swap<T&, U&> >::type::result_type
swap(T& t, U const& u)
{
   return Swap<T&, U&>()(t, const_cast<T&>(u));
}

template <typename T, typename U>
typename boost::enable_if<boost::mpl::and_<is_mutable_proxy<T>, is_mutable_proxy<U> >,
                          Swap<T&, U&> >::type::result_type
swap(T const& t, U const& u)
{
   return Swap<T&, U&>()(const_cast<T&>(t), const_cast<T&>(u));
}

// Transform
template <typename T, typename F, typename Enable = void>
struct Transform;  // primary version defined later, to handle identity functors

template <typename T, typename F, typename TInterface = typename interface<T>::type,
          typename Enable = void>
struct TransformInterface {};

#if 0
template <typename T, typename F>
inline
typename Transform<T, F>::result_type
transform(T const& x, F const& f)
{
   return Transform<T, F>()(x, f);
}
#endif

template <typename T, typename F, typename TInterface = typename interface<T>::type,
          typename Enable = void>
struct TransformRef;  // primary version defined later, to handle identity functors

template <typename T, typename F>
inline
typename Transform<typename proxy_lvalue<T>::type, F>::result_type
transform(T const& x, F const& f)
{
   return Transform<typename proxy_lvalue<T>::type, F>()(proxy_lvalue_cast(x), f);
}

template <typename T, typename F>
inline
typename Transform<T&, F>::result_type
transform(T& x, F const& f)
{
   return Transform<T&, F>()(x, f);
}

// ForwardIfDefined: if F is defined, then forward to T.  Else empty.
template <typename F, typename T, typename Enable = void>
struct ForwardIfDefined
{
};

template <typename F, typename T>
struct ForwardIfDefined<F, T, typename boost::enable_if<is_defined<F> >::type>
   : T {};

// TransformIfDefined - forwards to Transform<T, F> only if F is defined.
template <typename T, typename F, typename Enable = void>
struct TransformIfDefined : ForwardIfDefined<F, Transform<T, F> > {};

// TransformRefIfDefined - forwards to TransformRef<T, F> only if F is defined.
template <typename T, typename F, typename Enable = void>
struct TransformRefIfDefined : ForwardIfDefined<F, TransformRef<T, F> > {};

// Transform by an identity operator is a NOP

template <typename T, typename F, typename TInterface = typename interface<T>::type,
          typename Enable = void>
struct TransformIfIdentity : TransformInterface<T,F> {};

template <typename T, typename F, typename Ti>
struct TransformIfIdentity<T, F, Ti, typename boost::enable_if<is_identity<F> >::type>
{
   typedef boost::mpl::true_ identity;

   typedef typename make_const_reference<T>::type argument_type;
   typedef argument_type first_argument_type;
   typedef F second_argument_type;
   typedef argument_type result_type;
   result_type operator()(argument_type x) const { return x; }
   result_type operator()(first_argument_type x, F) const { return x; }
};

template <typename T, typename F, typename Ti>
struct TransformIfIdentity<T&, F, Ti, typename boost::enable_if<is_identity<F> >::type>
{
   typedef boost::mpl::true_ identity;

   typedef typename make_reference<T>::type argument_type;
   typedef argument_type first_argument_type;
   typedef F second_argument_type;
   typedef argument_type result_type;
   result_type operator()(argument_type x) const { return x; }
   result_type operator()(first_argument_type x, F) const { return x; }
};

template <typename T, typename F, typename Enable>
struct Transform : TransformIfIdentity<T, F> {};

//template <typename T, typename F, typename TI, typename Enable>
//struct TransformRef : TransformIfIdentity<T, F> {};

// unary_compose

template <typename T1, typename T2, typename Enable = void>
struct Compose
{
   typedef T1 first_argument_type;
   typedef T2 second_argument_type;
   typedef UnaryComposer<T1, T2> result_type;
   result_type operator()(T1 const& x, T2 const& y) const { return result_type(x,y); }
};

template <typename T1, typename T2>
inline
typename Compose<T1, T2>::result_type
compose(T1 const& x, T2 const& y)
{
   return Compose<T1, T2>()(x,y);
}

template <typename T>
struct Compose<T, T, typename boost::enable_if<is_involutary<T> >::type>
{
   typedef T first_argument_type;
   typedef T second_argument_type;
   typedef Identity<typename T::result_type> result_type;

   result_type operator()(T const&, T const&) const { return result_type(); }
};

template <typename T>
struct Compose<T, T, typename boost::enable_if<is_idempotent<T> >::type>
{
   typedef T first_argument_type;
   typedef T second_argument_type;
   typedef T result_type;

   result_type operator()(T const& x, T const&) const { return x; }
};

template <typename T1, typename T2>
struct Inverse<UnaryComposer<T1, T2>,
  typename boost::enable_if<
   boost::mpl::and_<exists<typename Inverse<T1>::result_type>,
                    exists<typename Inverse<T2>::result_type> > >::type>
  : public Compose<typename Inverse<T2>::result_type,
                    typename Inverse<T1>::result_type>
{
   typedef Compose<typename Inverse<T2>::result_type,
                   typename Inverse<T1>::result_type> base;

   Inverse(UnaryComposer<T1, T2> const& x) : base(inverse(x.second), inverse(x.first)) {}

};

// transform compose<f1, compose<f2, f3> > into compose<compose<f1, f2>, f3>.  We can
// always do this, as there is no associative rule here, both are exactly equivalent.

template <typename T1, typename T2, typename T3>
struct Compose<T1, UnaryComposer<T2, T3> >
{
   typedef Compose<T1, T2> FirstComposer;
   typedef typename arg_type<typename FirstComposer::result_type>::type T1T2Type;
   typedef Compose<T1T2Type, T3> SecondComposer;
   typedef typename SecondComposer::result_type result_type;
   typedef T1 first_argument_type;
   typedef UnaryComposer<T2, T3> second_argument_type;
   result_type operator()(first_argument_type const& x, second_argument_type const& y)
   {
      return compose(compose(x, y.first), y.second);
   }
};

// compose with identity functions

template <typename T, typename U>
struct Compose<T, U, typename boost::enable_if<
   boost::mpl::and_<
      is_identity<T>
    , boost::mpl::not_<is_identity<U> >
   >
>::type>
{
   typedef T first_argument_type;
   typedef U second_argument_type;
   typedef U result_type;
   U operator()(T const&, U const& y) const { return y; }
};

template <typename T, typename U>
struct Compose<T, U, typename boost::enable_if<
   boost::mpl::and_<
      boost::mpl::not_<is_identity<T> >
    , is_identity<U>
   >
>::type>
{
   typedef T first_argument_type;
   typedef U second_argument_type;
   typedef T result_type;
   T operator()(T const& x, U const&) const { return x; }
};

template <typename T, typename U>
struct Compose<T, U, typename boost::enable_if<
   boost::mpl::and_<
      is_identity<T>
    , is_identity<U>
   >
>::type>
{
   typedef T first_argument_type;
   typedef U second_argument_type;
   typedef Identity<typename T::result_type> result_type;
   result_type operator()(T const&, U const&) const { return result_type(); }
};

// compose<real, conj>

template <typename T>
struct Compose<Real<T>, Conj<T> >
{
   typedef Real<T> first_argument_type;
   typedef Conj<T> second_argument_type;
   typedef Real<T> result_type;
   result_type operator()(first_argument_type const& x, second_argument_type const&) const
   { return x; }
};

// compose<imag, conj>

template <typename T>
struct Compose<Imag<T>, Conj<T> >
{
   typedef Imag<T> first_argument_type;
   typedef Conj<T> second_argument_type;
   typedef Compose<Negate<typename arg_type<typename Imag<T>::result_type>::type>, Imag<T> > Composer;
   typedef typename Composer::result_type result_type;
   result_type operator()(first_argument_type const&, second_argument_type const&) const
   { return result_type(); }
};

//
// binary functions
//

// BinaryTransform, strictly only a binary function if the functor has
// a default ctor.

template <typename S, typename T, typename Func,
          typename SInterface = typename interface<S>::type,
          typename TInterface = typename interface<T>::type,
          typename Enable = void>
struct BinaryTransform
{
};

template <typename S, typename T, typename Func>
typename BinaryTransform<S, T, Func>::result_type
transform(S const& x, T const& y, Func const& f)
{
   return BinaryTransform<S, T, Func>()(x,y,f);
}

// EqualTo

template <typename S, typename T, typename SInterface = typename interface<S>::type,
          typename TInterface = typename interface<T>::type,
          typename Enable = void>
struct EqualToInterface
{
};

template <typename S, typename T, typename Enable = void>
struct EqualTo : EqualToInterface<S, T> {};

struct EqualToF
{
   template <typename S, typename T>
   struct result
   {
      typedef typename EqualTo<S, T>::result_type type;
   };

   template <typename S, typename T>
   typename EqualTo<S, T>::result_type operator()(S const& x, T const& y) const
   {
      return EqualTo<S, T>()(x,y);
   }
};

template <typename S, typename T>
inline
typename EqualTo<S, T>::result_type
equal_to(S const& x, T const& y)
{
   return EqualTo<S, T>()(x,y);
}

template <typename S, typename T>
inline
typename boost::disable_if<is_builtin<EqualTo<S, T> >,
                           typename EqualTo<S, T>::result_type>::type
operator==(S const& x, T const& y)
{
   return EqualTo<S, T>()(x,y);
}

template <typename S, typename T>
inline
typename boost::disable_if<is_builtin<EqualTo<S, T> >,
                           typename EqualTo<S, T>::result_type>::type
operator!=(S const& x, T const& y)
{
   return !EqualTo<S, T>()(x,y);
}

// Addition

template <typename S, typename T, typename SInterface = typename interface<S>::type,
          typename TInterface = typename interface<T>::type,
          typename Enable = void>
struct AdditionInterface
{
};

template <typename S, typename T, typename Enable = void>
struct Addition : AdditionInterface<S, T> {};

struct AdditionF
{
   template <typename S, typename T>
   struct result
   {
      typedef typename Addition<S, T>::result_type type;
   };

   template <typename S, typename T>
   typename Addition<S, T>::result_type operator()(S const& x, T const& y) const
   {
      return Addition<S, T>()(x,y);
   }
};

template <typename S, typename T>
inline
typename boost::disable_if<is_builtin<Addition<S, T> >,
                           typename Addition<S, T>::result_type>::type
operator+(S const& x, T const& y)
{
   return Addition<S, T>()(x,y);
}

// Subtraction

template <typename S, typename T, typename SInterface = typename interface<S>::type,
          typename TInterface = typename interface<T>::type,
          typename Enable = void>
struct SubtractionInterface {};

template <typename S, typename T, typename Enable = void>
struct Subtraction : SubtractionInterface<S, T> {};

struct SubtractionF
{
   template <typename S, typename T>
   struct result
   {
      typedef typename Subtraction<S, T>::result_type type;
   };

   template <typename S, typename T>
   typename Subtraction<S, T>::result_type operator()(S const& x, T const& y) const
   {
      return Subtraction<S, T>()(x,y);
   }
};

template <typename S, typename T>
inline
typename boost::disable_if<is_builtin<Subtraction<S, T> >,
                           typename Subtraction<S, T>::result_type>::type
operator-(S const& x, T const& y)
{
   return Subtraction<S, T>()(x,y);
}

// Multiplication

template <typename S, typename T, typename SInterface = typename interface<S>::type,
          typename TInterface = typename interface<T>::type, typename Enable = void>
struct MultiplicationInterface
{
};

template <typename S, typename T, typename Enable = void>
struct Multiplication : MultiplicationInterface<S, T> {};

struct MultiplicationF
{
   template <typename S, typename T>
   struct result
   {
      typedef typename Multiplication<S, T>::result_type type;
   };

   template <typename S, typename T>
   typename Multiplication<S, T>::result_type operator()(S const& x, T const& y) const
   {
      return Multiplication<S, T>()(x,y);
   }
};

template <typename S, typename T>
inline
typename boost::disable_if<is_builtin<Multiplication<S, T> >,
                           Multiplication<S, T> >::type::result_type
operator*(S const& x, T const& y)
{
   return Multiplication<S, T>()(x,y);
}

//
// direct_sum
//

template <typename S, typename T,
          typename SInterface = typename interface<S>::type,
          typename TInterface = typename interface<T>::type,
          typename Enable = void>
struct DirectSumInterface {};

template <typename S, typename T, typename Enable = void>
struct DirectSum : DirectSumInterface<S, T> {};

template <typename S, typename T>
inline
typename DirectSum<S, T>::result_type
direct_sum(S const& x, T const& y)
{
   return DirectSum<S, T>()(x,y);
}

template <typename S, typename T, typename U>
inline
typename DirectSum<typename DirectSum<S, T>::result_type,
                   U>::result_type
direct_sum(S const& x, T const& y, U const& z)
{
   return direct_sum(direct_sum(x,y), z);
}

//
// direct_product
//

template <typename S, typename T,
          typename SInterface = typename interface<S>::type,
          typename TInterface = typename interface<T>::type,
          typename Enable = void>
struct DirectProductInterface {};

template <typename S, typename T, typename Enable = void>
struct DirectProduct : DirectProductInterface<S, T> {};

template <typename S, typename T>
inline
typename DirectProduct<S, T>::result_type
direct_product(S const& x, T const& y)
{
   return DirectProduct<S, T>()(x,y);
}

//
// inner_prod
//

template <typename S, typename T, typename SInterface = typename interface<S>::type,
          typename TInterface = typename interface<T>::type,
          typename Enable = void>
struct InnerProdInterface
{
};

template <typename S, typename T, typename Enable = void>
struct InnerProd : InnerProdInterface<S, T> {};

struct InnerProdF
{
   template <typename S, typename T>
   struct result
   {
      typedef typename InnerProd<S, T>::result_type type;
   };

   template <typename S, typename T>
   typename InnerProd<S, T>::result_type operator()(S const& x, T const& y) const
   {
      return InnerProd<S, T>()(x,y);
   }
};

template <typename S, typename T>
inline
typename boost::disable_if<is_builtin<InnerProd<S, T> >,
                           typename InnerProd<S, T>::result_type>::type
inner_prod(S const& x, T const& y)
{
   return InnerProd<S, T>()(x,y);
}

//
// BindFirst and BindSecond
//
// Make a unary function out of a binary function, by binding a constant
// to one of the arguments.
//
// TODO: optimize for the case where the functor is empty
//

template <typename BinaryOp>
struct BinderFirst
{
   typedef BinaryOp binary_functor_type;
   typedef typename BinaryOp::first_argument_type   bound_type;
   typedef typename BinaryOp::second_argument_type  arg_type;
   //   typedef typename BinaryOp::value_type            value_type;
   typedef typename BinaryOp::result_type           result_type;

   typedef typename boost::add_reference<typename boost::add_const<bound_type>::type>::type
   bound_arg_type;

   typedef typename boost::add_reference<typename boost::add_const<arg_type>::type>::type
   argument_type;

   BinderFirst(bound_arg_type b_, binary_functor_type const& f_ = binary_functor_type())
     : b(b_), f(f_) {}

   result_type operator()(argument_type x) const { return f(b, x); }

   bound_type b;
   binary_functor_type f;
};

template <typename BinaryOp>
struct BinderSecond
{
   typedef BinaryOp binary_functor_type;
   typedef typename BinaryOp::first_argument_type arg_type;
   typedef typename BinaryOp::second_argument_type bound_type;
   //   typedef typename BinaryOp::value_type value_type;
   typedef typename BinaryOp::result_type result_type;

   typedef typename boost::add_reference<typename boost::add_const<bound_type>::type>::type
   bound_arg_type;

   typedef typename boost::add_reference<typename boost::add_const<arg_type>::type>::type
   argument_type;

   BinderSecond(bound_arg_type b_, binary_functor_type const& f_ = binary_functor_type())
     : b(b_), f(f_) {}

   BinderSecond(BinderSecond const& x) : b(x.b), f(x.f) {}

   result_type operator()(argument_type x) const { return f(x, b); }

   bound_type b;
   binary_functor_type f;
};

//
// Commute
//

template <typename F>
struct Commuter
{
   typedef typename F::second_argument_type first_argument_type;
   typedef typename F::first_argument_type second_argument_type;
   typedef typename F::result_type result_type;

   typedef is_stateless<F> stateless;

   Commuter() : f_() {}
   Commuter(F const& f) : f_(f) {}

   result_type operator()(first_argument_type x, second_argument_type y) const
   { return f_(y,x); }

   F f_;
};

template <typename F, typename Enable = void>
struct CommuteDefault
{
   typedef F argument_type;
   typedef Commuter<F> result_type;
   result_type operator()() const { return Commuter<F>(); }
   result_type operator()(F const& f) const { return Commuter<F>(f); }
};

template <typename F>
struct CommuteDefault<
   F
 , typename boost::enable_if<
      boost::mpl::and_<
         exists<typename F::commuted_function>
       , is_stateless<F>
      >
   >::type
>
{
   typedef F argument_type;
   typedef typename F::commuted_function result_type;
   result_type operator()() const { return result_type(); }
   result_type operator()(F const& f) const { return result_type(); }
};

template <typename F>
struct CommuteDefault<F, typename boost::enable_if<
   boost::mpl::and_<
      exists<typename F::commuted_function>
    , boost::mpl::not_<is_stateless<F> >
   >
>::type>
{
   typedef F argument_type;
   typedef typename F::commuted_function result_type;
   result_type operator()(F const& f) const { return result_type(f); }
};

template <typename F>
struct Commute : CommuteDefault<F> {};

template <typename F>
struct Commute<Commuter<F> >
{
   typedef Commuter<F> const& argument_type;
   typedef F result_type;
   F operator()(argument_type x) const { return x.f_; }
};

template <typename F>
inline
typename Commute<F>::result_type
commute(F const& f)
{
   return Commute<F>()(f);
}

//
// bind_first, bind_second free functions
//

template <typename T>
struct BindSecond
{
   typedef T const& first_argument_type;
   typedef typename T::second_argument_type second_argument_type;
   typedef BinderSecond<T> result_type;
   result_type operator()(first_argument_type f, second_argument_type y) const
   { return result_type(y,f); }
};

#if 1
template <typename T>
struct BindFirst
{
   typedef T const& first_argument_type;
   typedef typename T::first_argument_type second_argument_type;
   typedef BinderFirst<T> result_type;
   result_type operator()(first_argument_type f, second_argument_type y) const
   { return result_type(y,f); }
};
#else
template <typename T>
struct BindFirst
{
   typedef T const& first_argument_type;
   typedef typename T::first_argument_type second_argument_type;
   typedef typename BindSecond<typename Commute<T>::result_type>::result_type result_type;
   result_type operator()(first_argument_type f, second_argument_type y) const
   { return bind_second(commute(f), y); }
};
#endif

template <typename T, typename U>
inline
typename BindFirst<T>::result_type
bind_first(T const& f, U const& x)
{
   return BindFirst<T>()(f,x);
}

template <typename T, typename U>
inline
typename BindSecond<T>::result_type
bind_second(T const& f, U const& x)
{
   return BindSecond<T>()(f,x);
}

// equal

template <typename T, typename U, typename Tol = double,
          typename TInterface = typename interface<T>::type,
          typename UInterface = typename interface<U>::type>
struct EqualInterface {};

template <typename T, typename U, typename Tol = double>
struct Equal : EqualInterface<T, U, Tol>
{
   Equal() : EqualInterface<T, U, Tol>() {};
   Equal(Tol x) : EqualInterface<T, U, Tol>(x) {};
};

template <typename T, typename U, typename Tol>
inline
typename Equal<T,U>::result_type
equal(T const& x, U const& y, Tol t)
{
   return Equal<T,U,Tol>(t)(x,y);
}

template <typename T, typename U>
inline
typename Equal<T,U>::result_type
 equal(T const& x, U const& y)
{
   return Equal<T,U>()(x,y);
}

// stream insertion.
// The stream insertion operator << uses SFINAE on StreamInsert<T>::result_type.
// To enable streaming of a type, specialize StreamInsert, define
// result_type to be std::ostream&, and define operator() to do the
// actual streaming.  StreamInsert is a binary function.

template <typename T, typename TInterface = typename interface<T>::type>
struct StreamInsert
{
};

template <typename T>
inline
typename StreamInsert<T>::result_type
operator<<(std::ostream& out, T const& x)
{
   return StreamInsert<T>()(out, x);
}

//
// EvalExpression - evaluates an expression.
//

template <typename T, typename TInterface = typename interface<T>::type>
struct EvalExpression : Identity<T> {};

template <typename T>
inline
typename EvalExpression<T>::result_type
eval_expression(T const& x)
{
   return EvalExpression<T>()(x);
}

template <typename T>
inline
typename EvalExpression<T>::result_type
eval(T const& x)
{
   return EvalExpression<T>()(x);
}

// ZeroAll - zeros all elements.

template <typename T, typename TInterface = typename interface<T>::type,
          typename Enable = void>
struct ZeroAllInterface {};

template <typename T, typename Enable = void>
struct ZeroAll {};

template <typename T, typename Enable>
struct ZeroAll<T&, Enable> : ZeroAllInterface<T&, typename interface<T>::type> {};

template <typename T>
inline
typename ZeroAll<T&>::result_type
zero_all(T& v)
{
   return ZeroAll<T&>()(v);
}

template <typename T>
inline
typename boost::enable_if<is_mutable_proxy<T>, ZeroAll<T&>>::type::result_type
zero_all(T const& v)
{
   return ZeroAll<T&>()(const_cast<T&>(v));
}

// is_zero - function returns true if argument is equal to zero

template <typename T, typename TInterface = typename interface<T>::type>
struct IsZeroInterface
{
   typedef bool result_type;
   bool operator()(T const& x) const
   {
      return false;
   }
};

template <typename T, typename Enable> // primary is in interface.h
struct IsZero : IsZeroInterface<T> {};

template <typename T>
inline
typename IsZero<T>::result_type
is_zero(T const& x)
{
   return IsZero<T>()(x);
}

//
// assign.
//
// By default, Assign forwards to AssignInterface; this gives us two
// bites at the cherry to specialize on expression types; firstly
// Assign itself can be specialized on specific LHS and RHS types.
// The default Assign forwards to AssignInterface which allows us to
// specialize based on the interface.
// We probably only want to specialize based on the right hand side
// of Assign; specializing on LHS too would most likely lead to
// partial ordering complications.
//

template <typename LHS, typename RHS,
          typename LHSInterface = typename interface<LHS>::type,
          typename RHInterface = typename interface<RHS>::type>
struct AssignInterface {};

template <typename LHS, typename RHS>
struct Assign : public AssignInterface<LHS, RHS> { };

template <typename LHS, typename RHS>
inline
typename Assign<LHS&, RHS>::result_type
assign(LHS& x, RHS const& y)
{
   return Assign<LHS&, RHS>()(x,y);
}

// yet another experiment in how to handle proxies
template <typename LHS, typename RHS>
inline
typename boost::enable_if<is_mutable_proxy<LHS>,
                          Assign<LHS&, RHS> >::type::result_type
assign(LHS const& x, RHS const& y)
{
   return Assign<LHS&, RHS>()(const_cast<LHS&>(x),y);
}

template <typename LHS, typename RHS>
struct Assign<value_with_zero<LHS>&, RHS>
{
   typedef value_with_zero<LHS>& result_type;
   typedef value_with_zero<LHS>& first_argument_type;
   typedef RHS const& second_argument_type;

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      x = y;
      return x;
   }
};

template <typename LHS, typename RHS>
struct Assign<LHS&, value_with_zero<RHS> >
{
   typedef typename Assign<LHS&, RHS>::result_type result_type;
   typedef LHS& first_argument_type;
   typedef value_with_zero<RHS> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      if (y.is_zero())
         assign(x, zero_or_die<LHS>());
      else
         assign(x, y.get());
      return x;
   }
};

template <typename LHS, typename RHS>
struct Assign<value_with_zero<LHS>&, value_with_zero<RHS> >
{
   typedef value_with_zero<LHS>& result_type;
   typedef value_with_zero<LHS>& first_argument_type;
   typedef value_with_zero<RHS> const& second_argument_type;

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      if (!y.is_zero())
      {
         if (x.is_zero())
            assign(x, y.get());
         else
            assign(x.get(), y.get());
      }
      return x;
   }
};

// assign_copy

template <typename LHS, typename RHS, typename Enable = void>
struct AssignCopy {};

template <typename LHS, typename RHS>
struct AssignCopy<
   LHS&,
   RHS,
   typename boost::enable_if<
      boost::mpl::and_<
         boost::mpl::not_<is_noalias<RHS> >
       , boost::mpl::not_<is_noalias<LHS> >
       , is_defined<Assign<LHS&, typename make_value<RHS>::type> >
      >
   >::type
>
{
   typedef LHS& first_argument_type;
   typedef RHS const& second_argument_type;
   typedef typename Assign<LHS&, typename make_value<RHS>::type>::result_type result_type;

   result_type operator()(LHS& x, RHS const& y) const
   {
      typename make_value<RHS>::type Temp(y);
      return assign(x,Temp);
   }
};

template <typename LHS, typename RHS>
struct AssignCopy<
   LHS&,
   RHS,
   typename boost::enable_if<
      boost::mpl::and_<
         boost::mpl::or_<is_noalias<RHS>, is_noalias<LHS> >
       , is_defined<Assign<LHS&, RHS> >
      >
   >::type
> : Assign<LHS&, RHS> {};

template <typename LHS, typename RHS>
inline
typename AssignCopy<LHS&, RHS>::result_type
assign_copy(LHS& x, RHS const& y)
{
   return AssignCopy<LHS&, RHS>()(x,y);
}

// yet another experiment in how to handle proxies
template <typename LHS, typename RHS>
inline
typename boost::enable_if<is_mutable_proxy<LHS>,
                          AssignCopy<LHS&, RHS> >::type::result_type
assign_copy(LHS const& x, RHS const& y)
{
   return AssignCopy<LHS&, RHS>()(const_cast<LHS&>(x),y);
}

// add.

template <typename LHS, typename RHS,
          typename LHSInterface = typename interface<LHS>::type,
          typename RHInterface = typename interface<RHS>::type>
struct AddInterface
{
};

template <typename LHS, typename RHS>
struct Add : public AddInterface<LHS, RHS> { };

template <typename LHS, typename RHS>
inline
typename Add<LHS&, RHS>::result_type
add(LHS& x, RHS const& y)
{
   return Add<LHS&, RHS>()(x,y);
}

template <typename LHS, typename RHS>
inline
typename boost::enable_if<is_mutable_proxy<LHS>,
                          Add<LHS&, RHS> >::type::result_type
add(LHS const& x, RHS const& y)
{
   return Add<LHS&, RHS>()(const_cast<LHS&>(x),y);
}


// value_with_zero specializations

template <typename LHS, typename RHS>
struct Add<value_with_zero<LHS>&, RHS>
{
   typedef value_with_zero<LHS>& result_type;
   typedef value_with_zero<LHS>& first_argument_type;
   typedef RHS const& second_argument_type;

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      if (x.is_zero())
         assign(x, y);
      else
         add(x.get(), y);
      return x;
   }
};

template <typename LHS, typename RHS>
struct Add<LHS&, value_with_zero<RHS> >
{
   typedef typename Add<LHS&, RHS>::result_type result_type;
   typedef LHS& first_argument_type;
   typedef value_with_zero<RHS> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      if (!y.is_zero())
         add(x, y.get());
   }
};

template <typename LHS, typename RHS>
struct Add<value_with_zero<LHS>&, value_with_zero<RHS> >
{
   typedef void result_type;
   typedef value_with_zero<LHS>& first_argument_type;
   typedef value_with_zero<RHS> const& second_argument_type;

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      if (!y.is_zero())
      {
         if (x.is_zero())
            assign(x, y.get());
         else
            add(x.get(), y.get());
      }
   }
};

// add_copy

template <typename LHS, typename RHS, typename Enable = void>
struct AddCopyBase1 {};

template <typename LHS, typename RHS, typename Enable = void>
struct AddCopyBase2 : AddCopyBase1<LHS, RHS> {};

template <typename LHS, typename RHS, typename Enable = void>
struct AddCopy : AddCopyBase2<LHS, RHS> {};

template <typename LHS, typename RHS>
struct AddCopyBase1<
   LHS&,
   RHS,
   typename boost::enable_if<
      boost::mpl::and_<
         boost::mpl::not_<is_noalias<RHS> >
       , boost::mpl::not_<is_noalias<LHS> >
       , is_defined<Add<LHS&, typename make_value<RHS>::type> >
      >
   >::type
>
{
   typedef typename Add<LHS&, typename make_value<RHS>::type>::result_type result_type;
   typedef LHS& first_argument_type;
   typedef RHS const& second_argument_type;

   result_type operator()(LHS& x, RHS const& y) const
   {
      typename make_value<RHS>::type Temp(y);
      return add(x,Temp);
   }
};

// this is a bit of a hack because sometimes the AddCopy function
// returns a reference to the result (similar to operator+= convention)
// and sometimes it returns void.
template <typename T>
struct fudge
{
   static T const& apply(T const& x) { return x; }
};

template <>
struct fudge<void>
{
   template <typename T>
   static void apply(T const& x) { }
};

template <typename LHS, typename RHS>
struct AddCopyBase2<LHS&, value_with_zero<RHS> >
{
   typedef typename Add<LHS&, typename make_value<RHS>::type>::result_type result_type;
   typedef LHS& first_argument_type;
   typedef value_with_zero<RHS> const& second_argument_type;

   result_type
   operator()(LHS& x, value_with_zero<RHS> const& y) const
   {
      if (!y.is_zero()) add_copy(x, y.get());
      return fudge<result_type>::apply(x);
   }
};

template <typename LHS, typename RHS>
struct AddCopyBase1<
   LHS&,
   RHS,
   typename boost::enable_if<
      boost::mpl::and_<
         boost::mpl::or_<is_noalias<RHS>, is_noalias<LHS> >
       , is_defined<Add<LHS&, RHS> >
      >
   >::type
> : Add<LHS&, RHS> {};

template <typename LHS, typename RHS>
inline
typename AddCopy<LHS&, RHS>::result_type
add_copy(LHS& x, RHS const& y)
{
   return AddCopy<LHS&, RHS>()(x,y);
}

// yet another experiment in how to handle proxies
template <typename LHS, typename RHS>
inline
typename boost::enable_if<is_mutable_proxy<LHS>,
                          AddCopy<LHS&, RHS> >::type::result_type
add_copy(LHS const& x, RHS const& y)
{
   return AddCopy<LHS&, RHS>()(const_cast<LHS&>(x),y);
}

// subtract.
//

template <typename LHS, typename RHS,
          typename LHSInterface = typename interface<LHS>::type,
          typename RHInterface = typename interface<RHS>::type>
struct SubtractInterface
{
};

template <typename LHS, typename RHS>
struct Subtract : public SubtractInterface<LHS, RHS> { };

template <typename LHS, typename RHS>
inline
typename Subtract<LHS&, RHS>::result_type
subtract(LHS& x, RHS const& y)
{
   return Subtract<LHS&, RHS>()(x,y);
}

template <typename LHS, typename RHS>
typename boost::enable_if<is_mutable_proxy<LHS>,
                          Subtract<LHS&, RHS> >::type::result_type
subtract(LHS const& x, RHS const& y)
{
   return Subtract<LHS&, RHS>()(const_cast<LHS&>(x),y);
}

// subtract_copy

template <typename LHS, typename RHS, typename Enable = void>
struct SubtractCopyBase1 {};

template <typename LHS, typename RHS, typename Enable = void>
struct SubtractCopyBase2 : SubtractCopyBase1<LHS, RHS> {};

template <typename LHS, typename RHS, typename Enable = void>
struct SubtractCopy : SubtractCopyBase2<LHS, RHS> {};

template <typename LHS, typename RHS>
struct SubtractCopyBase1<
   LHS&,
   RHS,
   typename boost::enable_if<
      boost::mpl::and_<
         boost::mpl::not_<is_noalias<RHS> >
       , boost::mpl::not_<is_noalias<LHS> >
       , is_defined<Subtract<LHS&, typename make_value<RHS>::type> >
      >
   >::type
>
{
   typedef LHS& first_argument_type;
   typedef RHS const& second_argument_type;
   typedef typename Subtract<LHS&, typename make_value<RHS>::type>::result_type
   result_type;
   result_type operator()(LHS& x, RHS const& y) const
   {
      typename make_value<RHS>::type Temp(y);
      return subtract(x,Temp);
   }
};

template <typename LHS, typename RHS>
struct SubtractCopyBase2<LHS&, value_with_zero<RHS> >
{
   typedef typename Subtract<LHS&, typename make_value<RHS>::type>::result_type result_type;
   typedef LHS& first_argument_type;
   typedef value_with_zero<RHS> const& second_argument_type;

   result_type operator()(LHS& x, value_with_zero<RHS> const& y) const
   {
      if (!y.is_zero()) subtract_copy(x, y.get());
      return x;
   }
};

template <typename LHS, typename RHS>
struct SubtractCopyBase1<
   LHS&,
   RHS,
   typename boost::enable_if<
      boost::mpl::and_<
         boost::mpl::or_<is_noalias<RHS>, is_noalias<LHS> >
       , is_defined<Subtract<LHS&, RHS> >
      >
   >::type
> : Subtract<LHS&, RHS> {};

template <typename LHS, typename RHS>
inline
typename SubtractCopy<LHS&, RHS>::result_type
subtract_copy(LHS& x, RHS const& y)
{
   return SubtractCopy<LHS&, RHS>()(x,y);
}

// yet another experiment in how to handle proxies
template <typename LHS, typename RHS>
inline
typename boost::enable_if<is_mutable_proxy<LHS>,
                          SubtractCopy<LHS&, RHS> >::type::result_type
subtract_copy(LHS const& x, RHS const& y)
{
   return SubtractCopy<LHS&, RHS>()(const_cast<LHS&>(x),y);
}

// Multiply

template <typename LHS, typename RHS>
struct MultiplyDefault {};

template <typename LHS, typename RHS>
struct MultiplyDefault<LHS&, RHS>
{
   typedef LHS& result_type;
   typedef LHS& first_argument_type;
   typedef RHS const& second_argument_type;

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      return x *= y;
   }
};

template <typename LHS, typename RHS,
          typename LHSInterface = typename interface<LHS>::type,
          typename RHInterface = typename interface<RHS>::type,
          typename Enable = void>
struct MultiplyInterface : MultiplyDefault<LHS, RHS> {};

template <typename LHS, typename RHS, typename Enable = void>
struct Multiply : public MultiplyInterface<LHS, RHS> { };

template <typename LHS, typename RHS>
typename Multiply<LHS&, RHS>::result_type
multiply(LHS& x, RHS const& y)
{
   return Multiply<LHS&, RHS>()(x,y);
}

template <typename LHS, typename RHS>
typename boost::enable_if<is_mutable_proxy<LHS>, Multiply<LHS&, RHS> >::type::result_type
multiply(LHS const& x, RHS const& y)
{
   return Multiply<LHS&, RHS>()(const_cast<LHS&>(x), y);
}

// default_tolerance() and set_default_tolerance()
// a trick to keep everything in the header file
namespace Private
{

inline
double& implement_default_tolerance()
{
   static double tol = 64 * std::numeric_limits<double>::epsilon();
   return tol;
}

} // namespace Private

inline
double default_tolerance()
{
   return Private::implement_default_tolerance();
}

inline
void set_default_tolerance(double x)
{
   Private::implement_default_tolerance() = x;
}

} // namespace LinearAlgebra

#endif
