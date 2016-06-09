// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/vectoroperationsbase.h
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
  Defines various operators and expressions for Vector

  Created 2004-04-11 Ian McCulloch
*/

#if !defined(VECTOROPERATIONS_H_JHH7F7H74Y78YR7YIUFGO8QQ)
#define VECTOROPERATIONS_H_JHH7F7H74Y78YR7YIUFGO8QQ

#include "operations.h"
#include "iteroperations.h"
#include "vectorinterface.h"
//#include "slice.h"
#include "value_with_zero.h"
#if defined(USE_PSTREAM)
#include "pstream/pstream.h"
#endif

#include "scalar.h"

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>

namespace LinearAlgebra
{

// EvalExpression

template <typename T, typename V, typename U>
struct EvalExpression<T, LOCAL_VECTOR(V, U) > : Identity<T> {};

template <typename T, typename V, typename U>
struct EvalExpression<T, VECTOR_EXPRESSION(V, U) >
{
   typedef typename make_value<T>::type result_type;
   typedef T argument_type;
   result_type operator()(T const& x) const { return result_type(x); }
};

// Size

template <typename T, typename Interface = typename interface<T>::type>
struct SizeInterface {};


template <typename T, typename Enable = void>
struct Size : SizeInterface<T> {};

template <typename T>
inline
typename Size<T>::result_type
size(T const& x)
{
   return Size<T>()(x);
}

// Resize

template <typename T>
struct Resize {};

template <typename T>
inline
typename Resize<T&>::result_type
resize(T& v, size_type n)
{
   return Resize<T&>()(v, n);
}

template <typename T>
struct is_resizable : is_defined<Resize<T&> > {};

// try_resize - attempts a resize.  If T is resizable, then
// forward directly to resize().  Otherwise assert that the
// new size is the same as the old size and do nothing.

template <typename T>
inline
typename Resize<T&>::result_type
try_resize(T& v, size_type n)
{
   return Resize<T&>()(v, n);
}

template <typename T>
inline
typename boost::disable_if<is_defined<Resize<T&> >, void>::type
try_resize(T& v, size_type n)
{
   PRECONDITION_EQUAL(size(v), n);
}

// most unary operations are easy, we
// just apply the operation to each element.
// We use TransformIfDefined, which makes an empty
// function (ie. undefined) unless is_defined<operation>
// Similar effect to using boost::enable_if.
template <typename T, typename S, typename U>
struct NegateInterface<T, VECTOR_EXPRESSION(S,U) > 
   : TransformIfDefined<T, Negate<S> > 
{
};

template <typename T, typename S, typename U>
struct RealInterface<T, VECTOR_EXPRESSION(S,U) > : TransformIfDefined<T, Real<S> > {};

template <typename T, typename S, typename U>
struct RealInterface<T&, VECTOR_EXPRESSION(S,U) > : TransformIfDefined<T&, Real<S&> > {};

template <typename T, typename S, typename U>
struct ImagInterface<T, VECTOR_EXPRESSION(S,U) > : TransformIfDefined<T, Imag<S> > {};

template <typename T, typename S, typename U>
struct ImagInterface<T&, VECTOR_EXPRESSION(S,U) > : TransformRefIfDefined<T, Imag<S&> > {};

template <typename T, typename S, typename U>
struct ConjInterface<T, VECTOR_EXPRESSION(S,U) > : TransformIfDefined<T, Conj<S> > {};

template <typename T, typename S, typename U>
struct HermInterface<T, VECTOR_EXPRESSION(S,U) > : TransformIfDefined<T, Herm<S> > {};

template <typename T, typename S, typename U>
struct AbsInterface<T, VECTOR_EXPRESSION(S,U) > : TransformIfDefined<T, Abs<S> > {};

template <typename T, typename S, typename U>
struct ArgInterface<T, VECTOR_EXPRESSION(S,U) > : TransformIfDefined<T, Arg<S> > {};

template <typename T, typename S, typename U>
struct SinInterface<T, VECTOR_EXPRESSION(S,U) > : TransformIfDefined<T, Sin<S> > {};

template <typename T, typename S, typename U>
struct CosInterface<T, VECTOR_EXPRESSION(S,U) > : TransformIfDefined<T, Cos<S> > {};

template <typename T, typename S, typename U>
struct ExpInterface<T, VECTOR_EXPRESSION(S,U) > : TransformIfDefined<T, Exp<S> > {};

template <typename T, typename S, typename U>
struct TransposeInterface<T, VECTOR_EXPRESSION(S,U) > : TransformIfDefined<T, Transpose<S> > {};

//
// The interface for sparse vectors (but also works for dense of course)
//

// GetVectorElement

template <typename T,
	  typename TInterface = typename interface<T>::type>
struct GetVectorElementInterface { };

template <typename T, typename Enable = void>
struct GetVectorElement : GetVectorElementInterface<T> { };


template <typename T>
inline
typename boost::disable_if<is_mutable_proxy<T>, GetVectorElement<T> >::type::result_type
get_element(T const& v, size_type n)
{
   return GetVectorElement<T>()(v,n);
}

template <typename T>
inline
typename boost::enable_if<is_mutable_proxy<T>, GetVectorElement<T&> >::type::result_type
get_element(T const& v, size_type n)
{
   return GetVectorElement<T&>()(const_cast<T&>(v),n);
}

template <typename T>
inline
typename GetVectorElement<T&>::result_type
get_element(T& v, size_type n)
{
   return GetVectorElement<T&>()(v,n);
}

// default implementation

template <typename T, typename Enable = void>
struct GetVectorElementDefault {};

template <typename T>
struct GetVectorElementDefault<
   T&
 , typename boost::enable_if<exists<typename T::reference> >::type
>
{
   typedef typename T::reference result_type;
   typedef T& first_argument_type;
   typedef size_type second_argument_type;
   
   result_type operator()(first_argument_type v, second_argument_type n) const
   {
      return v[n];
   }
};

template <typename T>
struct GetVectorElementDefault<
   T
 , typename boost::enable_if<exists<typename T::const_reference> >::type
>
{
   typedef typename T::const_reference result_type;
   typedef T const& first_argument_type;
   typedef size_type second_argument_type;
   
   result_type operator()(first_argument_type v, second_argument_type n) const
   {
      return v[n];
   }
};

template <typename T, typename Tv, typename Ti>
struct GetVectorElementInterface<T, LOCAL_VECTOR(Tv, Ti)>
   : GetVectorElementDefault<T> {};

// defaults for stride/contig vectors

template <typename T, typename Tv, typename Ti>
struct GetVectorElementInterface<T, STRIDE_VECTOR(Tv, Ti)>
{
   typedef Tv const& result_type;
   typedef T const& first_argument_type;
   typedef size_type second_argument_type;
   
   result_type operator()(first_argument_type v, second_argument_type n) const
   {
      DEBUG_PRECONDITION(n < size(v))(n)(size(v));
      return *(data(v) + n*stride(v));
   }
};

template <typename T, typename Tv, typename Ti>
struct GetVectorElementInterface<T&, STRIDE_VECTOR(Tv, Ti)>
{
   typedef Tv& result_type;
   typedef T& first_argument_type;
   typedef size_type second_argument_type;
   
   result_type operator()(first_argument_type v, second_argument_type n) const
   {
      DEBUG_PRECONDITION(n < size(v))(n)(size(v));
      return *(data(v) + n*stride(v));
   }
};

template <typename T, typename Tv, typename Ti>
struct GetVectorElementInterface<T, CONTIGUOUS_VECTOR(Tv, Ti)>
{
   typedef Tv const& result_type;
   typedef T const& first_argument_type;
   typedef size_type second_argument_type;
   
   result_type operator()(first_argument_type v, second_argument_type n) const
   {
      DEBUG_PRECONDITION(n < size(v))(n)(size(v));
      return *(data(v) + n);
   }
};

template <typename T, typename Tv, typename Ti>
struct GetVectorElementInterface<T&, CONTIGUOUS_VECTOR(Tv, Ti)>
{
   typedef Tv& result_type;
   typedef T& first_argument_type;
   typedef size_type second_argument_type;
   
   result_type operator()(first_argument_type v, second_argument_type n) const
   {
      DEBUG_PRECONDITION(n < size(v))(n)(size(v));
      return *(data(v) + n);
   }
};

// SetElement

template <typename T, typename Value,
	  typename TInterface = typename interface<T>::type, typename Enable = void>
struct SetElement { };

template <typename T, typename Value>
inline
typename SetElement<T&, Value>::result_type
set_element(T& v, size_type n, Value const& x)
{
   SetElement<T&, Value>()(v, n, x);
}

// FIXME: this should be a separate function, rather than delegate to SetElement
template <typename T, typename Value>
inline
typename SetElement<T&, Value>::result_type
set_new_element(T& v, size_type n, Value const& x)
{
   SetElement<T&, Value>()(v, n, x);
}

// AddElement

template <typename T, typename Value,
	  typename TInterface = typename interface<T>::type, typename Enable = void>
struct AddElement { };

template <typename T, typename Value>
inline
typename AddElement<T&, Value>::result_type
add_element(T& v, size_type n, Value const& x)
{
   AddElement<T&, Value>()(v, n, x);
}

// SubtractElement

template <typename T, typename Value,
	  typename TInterface = typename interface<T>::type, typename Enable = void>
struct SubtractElement { };

template <typename T, typename Value>
inline
typename SubtractElement<T&, Value>::result_type
subtract_element(T& v, size_type n, Value const& x)
{
   SubtractElement<T&, Value>()(v, n, x);
}

// ZeroElement

template <typename T, typename TInterface = typename interface<T>::type, 
	  typename Enable = void>
struct ZeroElement { };

template <typename T>
inline
typename ZeroElement<T&>::result_type
zero_element(T& v, size_type n)
{
   ZeroElement<T&>()(v, n);
}

// iterate_at

template <typename T, typename Ti = typename interface<T>::type>
struct VectorIterateAtInterface {};

template <typename T, typename Enable = void>
struct VectorIterateAt : VectorIterateAtInterface<T> {};

template <typename T>
inline
typename VectorIterateAt<T>::result_type
iterate_at(T const& m, size_type i)
{
   return VectorIterateAt<T>()(m, i);
}

template <typename T>
inline
typename VectorIterateAt<T&>::result_type
iterate_at(T& m, size_type i)
{
   return VectorIterateAt<T&>()(m, i);
}

// default implementation for dense
template <typename T, typename Tv, typename Ti>
struct VectorIterateAt<T, DENSE_VECTOR(Tv, Ti)>
{
   typedef typename const_iterator<T>::type result_type;
   typedef T const& first_argument_type;
   typedef size_type second_argument_type;

   result_type operator()(first_argument_type m, second_argument_type i) const
   {
      typename const_iterator<T>::type I = iterate(m);
      I += i;
      return I;
   }
};

// sum

template <typename T, typename TV, typename U>
struct Sum<T, VECTOR_EXPRESSION(TV, U)>
   : Sum<typename EvalExpression<T>::result_type> {};

template <typename T, typename TV, typename U>
struct Sum<T, LOCAL_VECTOR(TV, U)>
{
   typedef typename make_value_with_zero<TV>::type result_type;
   //   typedef TV result_type;
   typedef T argument_type;

   result_type operator()(T const& x) const
   {
      return iter_sum(iterate(x));
   }
};

// norm_1

template <typename T, typename TV, typename U>
struct Norm1<T, VECTOR_EXPRESSION(TV, U)>
   : Norm1<typename EvalExpression<T>::result_type> {};

template <typename T, typename TV, typename U>
struct Norm1<T, LOCAL_VECTOR(TV, U)>
{
   typedef typename Norm1<TV>::result_type result_type;
   typedef T argument_type;
   result_type operator()(T const& x) const
   {
      return iter_norm_1(iterate(x));
   }
};

// norm_2

template <typename T, typename TV, typename U>
struct Norm2Sq<T, VECTOR_EXPRESSION(TV, U)>
   : Norm2Sq<typename EvalExpression<T>::result_type> {};

template <typename T, typename TV, typename U>
struct Norm2Sq<T, LOCAL_VECTOR(TV, U)>
{
   typedef typename Norm2Sq<TV>::result_type result_type;
   typedef T argument_type;
   result_type operator()(T const& x) const
   {
      return iter_norm_2_sq(iterate(x));
   }
};

// norm_inf

template <typename T, typename TV, typename U>
struct NormInf<T, VECTOR_EXPRESSION(TV, U)>
   : NormInf<typename EvalExpression<T>::result_type> {};

template <typename T, typename TV, typename U>
struct NormInf<T, LOCAL_VECTOR(TV, U)>
{
   typedef typename NormInf<TV>::result_type result_type;
   typedef T argument_type;
   result_type operator()(T const& x) const
   {
      return iter_norm_inf(iterate(x));
   }
};

// norm_frob

template <typename T, typename TV, typename U>
struct NormFrobSq<T, VECTOR_EXPRESSION(TV, U)>
   : NormFrobSq<typename EvalExpression<T>::result_type> {};

template <typename T, typename TV, typename U>
struct NormFrobSq<T, LOCAL_VECTOR(TV, U)>
{
   typedef typename NormFrobSq<TV>::result_type result_type;
   typedef T argument_type;
   result_type operator()(T const& x) const
   {
      return iter_norm_frob_sq(iterate(x));
   }
};

// min
// TODO: add a comparison function to min & max

template <typename T, typename TInterface = typename interface<T>::type>
struct Min {};

template <typename T>
inline
typename Min<T>::result_type
min(T const& x)
{
   return Min<T>()(x);
}

// default implementation

template <typename T, typename Tv, typename Ti>
struct Min<T, INJECTIVE_VECTOR(Tv, Ti)>
{
   typedef typename GetVectorElement<T>::result_type result_type;
   typedef T const& argument_type;
   result_type operator()(argument_type v) const
   {
      return *iter_min(iterate(v));
   }
};

template <typename T, typename Tv, typename Ti>
struct Min<T, DENSE_VECTOR(Tv, Ti)>
{
   typedef typename GetVectorElement<T>::result_type result_type;
   typedef T const& argument_type;
   result_type operator()(argument_type v) const
   {
      return *iter_min(iterate(v));
   }
};

// max

template <typename T, typename TInterface = typename interface<T>::type>
struct Max {};

template <typename T>
inline
typename Max<T>::result_type
max(T const& x)
{
   return Max<T>()(x);
}

// default implementation

template <typename T, typename Tv, typename Ti>
struct Max<T, INJECTIVE_VECTOR(Tv, Ti)>
{
   typedef typename GetVectorElement<T>::result_type result_type;
   typedef T const& argument_type;
   result_type operator()(argument_type v) const
   {
      return *iter_max(iterate(v));
   }
};

template <typename T, typename Tv, typename Ti>
struct Max<T, DENSE_VECTOR(Tv, Ti)>
{
   typedef typename GetVectorElement<T>::result_type result_type;
   typedef T const& argument_type;
   result_type operator()(argument_type v) const
   {
      return *iter_max(iterate(v));
   }
};

// binary operators

// inner_prod

template <typename S, typename T, 
	  typename Nested = InnerProd<typename interface<S>::value_type, 
				      typename interface<T>::value_type>,
	  typename SInterface = typename interface<S>::type, 
	  typename TInterface = typename interface<T>::type,
	  typename Enable = void>
struct VectorInnerProd
{
};

// overload of inner_prod for 3-argument version
template <typename S, typename T, typename Func>
inline
typename VectorInnerProd<S, T, Func>::result_type
inner_prod(S const& x, T const& y, Func const& f)
{
   return VectorInnerProd<S, T, Func>(f)(x,y);
}

template <typename S, typename T, typename Sv, typename Si, typename Tv, typename Ti>
struct InnerProdInterface<S, T, ANY_VECTOR(Sv, Si), ANY_VECTOR(Tv, Ti)>
   : VectorInnerProd<S, T> {};

template <typename S, typename T, typename Func, 
	  typename Sv, typename Si, typename Tv, typename Ti>
struct VectorInnerProd<S, T, Func, VECTOR_EXPRESSION(Sv, Si), VECTOR_EXPRESSION(Tv, Ti)>
   : VectorInnerProd<typename EvalExpression<S>::result_type,
		     typename EvalExpression<T>::result_type,
		     Func> 
{
   typedef VectorInnerProd<typename EvalExpression<S>::result_type,
      typename EvalExpression<T>::result_type,
      Func>  base;
   VectorInnerProd() {}
   VectorInnerProd(Func const& f) : base(f) {}
};

template <typename S, typename T, typename Func, 
	  typename Sv, typename Si, typename Tv, typename Ti>
struct VectorInnerProd<S, T, Func, LOCAL_VECTOR(Sv, Si), LOCAL_VECTOR(Tv, Ti)>
{
   typedef typename make_value_with_zero<typename Func::result_type>::type result_type;
   typedef S const& first_argument_type;
   typedef T const& second_argument_type;

   VectorInnerProd() {}
   VectorInnerProd(Func const& f) : f_(f) {}

   typedef Func const& third_argument_type;

   result_type operator()(S const& x, T const& y) const
   { return iter_inner_prod(iterate(x), iterate(y), f_); }

   result_type operator()(S const& x, T const& y, Func const& f) const
   { return iter_inner_prod(iterate(x), iterate(y), f); }

   Func f_;
};

// parallel_prod

template <typename S, typename T, 
	  typename Nested = Multiplication<typename interface<S>::value_type, 
					   typename interface<T>::value_type>,
	  typename SInterface = typename interface<S>::type, 
	  typename TInterface = typename interface<T>::type,
	  typename Enable = void>
struct ParallelProdInterface {};

template <typename S, typename T, 
	  typename Nested = Multiplication<typename interface<S>::value_type, 
					   typename interface<T>::value_type>,
	  typename Enable = void>
struct ParallelProd : ParallelProdInterface<S, T, Nested> 
{
   ParallelProd() {}
   ParallelProd(Nested const& n) : ParallelProdInterface<S, T, Nested>(n) {}
};


template <typename S, typename T, typename Nested,
          typename Sv, typename Si,
          typename Tv, typename Ti>
struct ParallelProdInterface<S, T, Nested, ANY_VECTOR(Sv, Si), ANY_VECTOR(Tv, Ti)>
   : VectorInnerProd<S, T, Nested> 
{
   ParallelProdInterface() {}
   ParallelProdInterface(Nested const& n) : VectorInnerProd<S, T, Nested>(n) {}
};

template <typename S, typename T, typename Func>
inline
typename ParallelProd<S, T, Func>::result_type
parallel_prod(S const& x, T const& y, Func const& f)
{
   //   typedef typename boost::mpl::print<typename Func::result_type>::type dummy;
   return ParallelProd<S, T, Func>(f)(x,y);
}

template <typename S, typename T>
inline
typename ParallelProd<S, T>::result_type
parallel_prod(S const& x, T const& y)
{
   return ParallelProd<S, T>()(x,y);
}

struct ParallelProdF
{
   template <typename Args> struct result {};

   template <typename T1, typename T2>
   struct result<ParallelProdF(T1,T2)>
   {
      typedef typename ParallelProd<T1,T2>::result_type type;
   };

   template <typename T1, typename T2>
   typename ParallelProd<T1,T2>::result_type
   operator()(T1 const& x, T2 const& y) const
   {
      return ParallelProd<T1,T2>()(x,y);
   }
};

// equal_to

template <typename S, typename T, typename V1, typename U1, typename V2, typename U2>
struct EqualToInterface<S, T, DENSE_VECTOR(V1, U1), DENSE_VECTOR(V2, U2) >
{
   typedef bool result_type;
   typedef S first_argument_type;
   typedef T second_argument_type;

   result_type operator()(S const& x, T const& y) const
   {
      return iter_equal_to(iterate(x), iterate(y));
   }
};

template <typename S, typename T, typename V1, typename U1, typename V2, typename U2>
struct EqualToInterface<S, T, INJECTIVE_VECTOR(V1, U1), INJECTIVE_VECTOR(V2, U2) >
{
   typedef bool result_type;
   typedef S first_argument_type;
   typedef T second_argument_type;

   result_type operator()(S const& x, T const& y) const
   {
      return iter_equal_to(iterate(x), iterate(y));
   }
};

// addition

template <typename S, typename T, typename Sv, typename Si, typename Tv, typename Ti>
struct AdditionInterface<S, T, VECTOR_EXPRESSION(Sv, Si), VECTOR_EXPRESSION(Tv, Ti)>
   : BinaryTransform<S, T, Addition<Sv, Tv> > { };

// subtraction

template <typename S, typename T, typename Sv, typename Si, typename Tv, typename Ti>
struct SubtractionInterface<S, T, VECTOR_EXPRESSION(Sv, Si), VECTOR_EXPRESSION(Tv, Ti)>
   : BinaryTransform<S, T, Subtraction<Sv, Tv> > { };


#if 0
// old stuff, new version of equal doesn't yet se interface matching
template <typename S1, typename D1, typename S2, typename D2>
inline
bool equal(VectorConstExpression<S1, D1> const& x, VectorConstExpression<S2, D2> const& y, 
	   double tol = default_tolerance())
{
   return x.size() == y.size() && norm_inf(x-y) <= 
      tol + 2 * std::numeric_limits<double>::epsilon() * (norm_inf(x) + norm_inf(y));
}
#endif

// DirectSum

template <typename S, typename T, 
	  typename SInterface = typename interface<S>::type,
	  typename TInterface = typename interface<T>::type>
struct VectorDirectSum {};

template <typename S, typename T, typename Sv, typename Si, typename Tv, typename Ti>
struct DirectSumInterface<S, T, ANY_VECTOR(Sv, Si), ANY_VECTOR(Tv, Ti)>
   : VectorDirectSum<S, T> {};

//
// generate
//
// set elements of a vector according to a generating function
//

template <typename Iter, typename Func>
void
generate_iter(Iter i, Func f)
{
   while (i)
   {
      *i = f();
      ++i;
   }
}

template <typename T, typename Func>
void
generate(T& M, Func f = Func())
{
   generate_iter(iterate(M), f);
}

template <typename Func>
Vector<typename make_value<typename Func::result_type>::type>
generate(size_type Size, Func f = Func())
{
   Vector<typename make_value<typename Func::result_type>::type> Result(Size);
   generate(Result, f);
   return Result;
}

template <typename Scalar>
Vector<typename make_value<Scalar>::type>
generate(size_type Size, Scalar(&f)())
{
   Vector<typename make_value<Scalar>::type> Result(Size);
   generate(Result, f);
   return Result;
}

// NNZ - number of non-zero elements, ordinary function

template <typename T, typename TInterface = typename interface<T>::type>
struct NNZInterface { };

template <typename T,  typename Enable = void>
struct NNZ : NNZInterface<T> { };

template <typename T>
inline
typename NNZ<T>::result_type
nnz(T const& x)
{
   return NNZ<T>()(x);
}

//
// Fill
//

template <typename T, typename V, typename TInterface = typename interface<T>::type>
struct Fill
{
};

template <typename T, typename V>
typename Fill<T&, V>::result_type
fill(T& x, V const& y)
{
   return Fill<T&,V>()(x,y);
}

template <typename T, typename V>
typename boost::enable_if<is_proxy<T>, Fill<T&, V> >::type::result_type
fill(T const& x, V const& y)
{
   return Fill<T&,V>()(const_cast<T&>(x),y);
}

// default implementation

template <typename T, typename V, typename Tv, typename Ti>
struct Fill<T&, V, DENSE_VECTOR(Tv, Ti)>
{
   typedef void result_type;
   typedef T& first_argument_type;
   typedef V second_argument_type;
   void operator()(T& v, V const& x) const
   {
      iter_fill(iterate(v), x);
   }
};

// Swap

template <typename T, typename U, typename V, typename Ti, typename Ui>
struct SwapInterface<T&, U&, DENSE_VECTOR(V, Ti), DENSE_VECTOR(V, Ui)>
{
   typedef void result_type;
   typedef T& first_argument_type;
   typedef U& second_argument_type;
   result_type operator()(T& t, U& u) const
   {
      Vector<V> Temp(t);
      assign(t, u);
      assign(u, Temp);
   }
};

//
// square-bracket operator
//

template <typename T, typename Arg, 
	  typename Ti = typename interface<T>::type,
	  typename Argi = typename interface<Arg>::type>
struct VectorBracketInterface { };

template <typename T, typename Arg, typename Enable = void>
struct VectorBracket : VectorBracketInterface<T, Arg> {};

// 
template <typename T, typename Arg>
struct VectorBracket<
   T
 , Arg
 , typename boost::enable_if<
      boost::mpl::and_<
         boost::is_integral<Arg>
       , exists<typename GetVectorElement<T>::result_type>
       , boost::mpl::not_<is_builtin<GetVectorElement<T> > >
      >
   >::type
> : GetVectorElement<T> {};

//
// multiplication
//

// vector * scalar

template <typename S, typename T,
	  typename NestedMultiply = Multiplication<typename interface<S>::value_type, T>,
	  typename SInterface = typename interface<S>::type, 
	  typename Enable = void>
struct VectorScalarMultiplication
{
};

template <typename S, typename T, typename SV, typename U, typename Mult>
struct VectorScalarMultiplication<S, T, Mult, ANY_VECTOR(SV, U) >
//   : Transform<S, typename BindSecond<Mult>::result_type>
{
   typedef typename BindSecond<Mult>::result_type Func;
   typedef S first_argument_type;
   //   typedef typename Mult::second_argument_type second_argument_type;
   typedef T second_argument_type;
   typedef typename Transform<S, Func>::result_type result_type;

   result_type operator()(S const& v, T const& x) const
   {
      return transform(v, bind_second(Mult(), x));
   }
};

// we forward to VectorScalarMultiplication<S, typename interface<T>::value_type> 
// rather than using T directly, to handle the case where T is a ScalarProxy.

template <typename S, typename T, typename SV, typename SI>
struct MultiplicationInterface<S, T, ANY_VECTOR(SV, SI), AnyScalar<T> >
   : VectorScalarMultiplication<S, typename interface<T>::value_type> 
{
};

// scalar * vector

template <typename S, typename T,
	  typename NestedMultiply = Multiplication<S, typename interface<T>::value_type>,
	  typename TInterface = typename interface<T>::type, 
	  typename Enable = void>
struct ScalarVectorMultiplication
{
};

template <typename S, typename T, typename TV, typename U, typename Mult>
struct ScalarVectorMultiplication<S, T, Mult, ANY_VECTOR(TV, U) >
//   : Transform<T, typename BindFirst<Mult>::result_type>
{
   typedef typename BindFirst<Mult>::result_type Func;
   typedef S first_argument_type;
   typedef T second_argument_type;
   typedef typename Transform<T, Func>::result_type result_type;

   result_type operator()(S const& x, T const& v) const
   {
      return transform(v, bind_first(Mult(), x));
   }
};

template <typename S, typename T, typename TV, typename TI>
struct MultiplicationInterface<S, T, AnyScalar<S>, ANY_VECTOR(TV, TI)>
   : ScalarVectorMultiplication<typename interface<S>::value_type, T> 
{
};

// stream output

// sparse vectors

template <typename T, typename S, typename U>
struct StreamInsert<T, COMPRESSED_VECTOR(S, U) >
{
   typedef std::ostream& result_type;
   typedef std::ostream& first_argument_type;
   typedef T second_argument_type;
   result_type operator()(std::ostream& out, T const& x) const;
};

template <typename T, typename S, typename U>
std::ostream&
StreamInsert<T, COMPRESSED_VECTOR(S, U) >::operator()(std::ostream& out, T const& x) const
{
   typename const_iterator<T>::type I = iterate(x);
   if (!I)
   {
      return out << "(empty)";
   }
   // else
   out << "((index=" << I.index() << ", value=" << *I << ")";
   ++I;
   while (I)
   {
      out << ", (index=" << I.index() << ", value=" << *I << ")";
      ++I;
   }
   out << ')';
   return out;
}

// dense vectors

template <typename T, typename S, typename U>
struct StreamInsert<T, DENSE_VECTOR(S,U)>
{
   typedef std::ostream& result_type;
   typedef std::ostream& first_argument_type;
   typedef T second_argument_type;
   result_type operator()(std::ostream& out, T const& x) const;
};

template <typename T, typename S, typename U>
std::ostream&
StreamInsert<T, DENSE_VECTOR(S,U) >::operator()(std::ostream& out, T const& x) const
{
   typename const_iterator<T>::type I = iterate(x);

   //typedef typename boost::mpl::print< typename const_iterator<T>::type>::type dummy;

   if (!I)
   {
      return out << "(empty)";
   }
   // else
   out << '(' << *I;
   ++I;
   while (I)
   {
      out << ", " << *I;
      ++I;
   }
   out << ')';
   return out;
}

// expressions

// some kind of compose metafunction would be nice here
template <typename T, typename S, typename U>
struct StreamInsert<T, VECTOR_EXPRESSION(S,U)>
{
   typedef std::ostream& result_type;
   typedef std::ostream& first_argument_type;
   typedef T second_argument_type;
   result_type operator()(std::ostream& out, T const& x) const
   {
      return out << eval_expression(x);
   }
};

//
// Data, Stride
//
// For STRIDE_VECTOR types, we have some additional functions,
// data(x) returns a pointer to the array (with Ref and non-Ref variants),
// stride(x) returns the stride.
//

template <typename T, typename TInterface = typename interface<T>::type>
struct DataInterface
{
};

template <typename T>
struct Data : DataInterface<T> {};

template <typename T>
typename Data<T>::result_type
data(T const& x)
{
   return Data<T>()(x);
}

template <typename T>
typename Data<T&>::result_type
data(T& x)
{
   return Data<T&>()(x);
}

// default implementation

template <typename T, typename Tv, typename Ti>
struct DataInterface<T, STRIDE_VECTOR(Tv, Ti)>
{
   typedef Tv const* result_type;
   typedef T argument_type;
   result_type operator()(T const& x) const { return x.data(); }
};

template <typename T, typename Tv, typename Ti>
struct DataInterface<T&, STRIDE_VECTOR(Tv, Ti)>
{
   typedef Tv* result_type;
   typedef T& argument_type;
   result_type operator()(T& x) const { return x.data(); }
};

// stride

template <typename T, typename TInterface = typename interface<T>::type>
struct Stride
{
};

template <typename T>
typename Stride<T>::result_type
stride(T const& x)
{
   return Stride<T>()(x);
}

// default implementation

template <typename T, typename Tv, typename Ti>
struct Stride<T, STRIDE_VECTOR(Tv, Ti)>
{
   typedef difference_type result_type;
   typedef T argument_type;
   result_type operator()(T const& x) const { return x.stride(); }
};

template <typename T, typename Tv, typename Ti>
struct Stride<T, CONTIGUOUS_VECTOR(Tv, Ti)>
{
   typedef difference_type result_type;
   typedef T argument_type;
   result_type operator()(T const&) const { return 1; }
};

// equal

template <typename T, typename U, typename TolType, 
	  typename Tv, typename Ti,
	  typename Uv, typename Ui>
struct EqualInterface<T, U, TolType, ANY_VECTOR(Tv,Ti), ANY_VECTOR(Uv, Ui)>
{
   typedef bool result_type;
   typedef T const& first_argument_type;
   typedef U const& second_argument_type;

   EqualInterface(TolType Tol = default_tolerance()) : Tol_(Tol) {}

   bool operator()(T const& x, U const& y) const
   {
      return (size(x) == size(y)) && norm_inf(x-y) 
	 <= Tol_ + 2 * std::numeric_limits<TolType>::epsilon() * (norm_inf(x) + norm_inf(y));
   }

   private:
      TolType Tol_;
};

// assign

template <typename LHS, typename RHS>
struct AssignExpression {};

template <typename LHS, typename RHS, typename S1, typename U1, typename S2, typename U2>
struct AssignInterface<LHS, RHS, VECTOR_EXPRESSION(S1, U1), VECTOR_EXPRESSION(S2, U2) >
   : AssignExpression<LHS, RHS>
{
};

// make assignment of a non-expression RHS fail by default
template <typename LHS, typename RHS, typename S1, typename U1, typename S2, typename U2>
struct AssignInterface<LHS, RHS, VECTOR_EXPRESSION(S1, U1), LOCAL_VECTOR(S2, U2) >
{
};

template <typename LHS, typename RHS, typename S1, typename U1, typename S2, typename U2>
struct AssignInterface<LHS&, RHS, DENSE_VECTOR(S1, U1), DENSE_VECTOR(S2, U2) >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      try_resize(x, size(y));
      DEBUG_PRECONDITION_EQUAL(size(x), size(y));
      iter_assign(Iterate<LHS&>()(x), Iterate<RHS>()(y));
   }
};

template <typename LHS, typename RHS, typename S1, typename U1, typename S2, typename U2>
struct AssignInterface<LHS&, RHS, DENSE_VECTOR(S1, U1), COMPRESSED_VECTOR(S2, U2) >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      try_resize(x, size(y));
      DEBUG_PRECONDITION_EQUAL(size(x), size(y));
      zero_all(x);
      typename Iterate<RHS>::result_type r = iterate(y);
      typename Iterate<LHS&>::result_type l = iterate(x);
      while (r)
      {
	 l[r.index()] += *r;
	 ++r;
      }
   }
};

template <typename LHS, typename RHS, typename S1, typename U1, typename S2, typename U2>
struct AssignInterface<LHS&, RHS, DENSE_VECTOR(S1, U1), INJECTIVE_VECTOR(S2, U2) >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      try_resize(x, size(y));
      DEBUG_PRECONDITION_EQUAL(size(x), size(y));
      zero_all(x);
      typename Iterate<RHS>::result_type r = iterate(y);
      typename Iterate<LHS&>::result_type l = iterate(x);
      while (r)
      {
	 l[r.index()] = *r;
	 ++r;
      }
   }
};

template <typename LHS, typename RHS, typename S1, typename U1, typename S2, typename U2>
struct AssignInterface<LHS&, RHS, COMPRESSED_VECTOR(S1, U1), COMPRESSED_VECTOR(S2, U2) >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      try_resize(x, size(y));
      DEBUG_PRECONDITION_EQUAL(size(x), size(y));
      zero_all(x);
      typename Iterate<RHS>::result_type r = iterate(y);
      while (r)
      {
	 add_element(x, r.index(), *r);
	 ++r;
      }
   }
};

template <typename LHS, typename RHS, typename S1, typename U1, typename S2, typename U2>
struct AssignInterface<LHS&, RHS, INJECTIVE_VECTOR(S1, U1), INJECTIVE_VECTOR(S2, U2) >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      try_resize(x, size(y));
      DEBUG_PRECONDITION_EQUAL(size(x), size(y));
      zero_all(x);
      typename Iterate<RHS>::result_type r = iterate(y);
      while (r)
      {
	 set_new_element(x, r.index(), *r);
	 ++r;
      }
   }
};

// add

template <typename LHS, typename RHS>
struct AddExpression
{
};

template <typename LHS, typename RHS, typename S1, typename U1, typename S2, typename U2>
struct AddInterface<LHS, RHS, VECTOR_EXPRESSION(S1, U1), VECTOR_EXPRESSION(S2, U2) >
   : AddExpression<LHS, RHS>
{
};

template <typename LHS, typename RHS, typename S1, typename U1, typename S2, typename U2>
struct AddInterface<LHS, RHS, VECTOR_EXPRESSION(S1, U1), LOCAL_VECTOR(S2, U2) >
{
};

template <typename LHS, typename RHS, typename S1, typename U1, typename S2, typename U2>
struct AddInterface<LHS&, RHS, DENSE_VECTOR(S1, U1), DENSE_VECTOR(S2, U2) >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      DEBUG_PRECONDITION_EQUAL(size(x), size(y));
      iter_add(iterate(x), iterate(y));
   }
};

template <typename LHS, typename RHS, typename S1, typename U1, typename S2, typename U2>
struct AddInterface<LHS&, RHS, DENSE_VECTOR(S1, U1), COMPRESSED_VECTOR(S2, U2) >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      DEBUG_PRECONDITION_EQUAL(size(x), size(y));
      typename Iterate<RHS>::result_type r = iterate(y);
      typename Iterate<LHS&>::result_type l = iterate(x);
      while (r)
      {
	 l[r.index()] += *r;
	 ++r;
      }
   }
};

template <typename LHS, typename RHS, typename S1, typename U1, typename S2, typename U2>
struct AddInterface<LHS&, RHS, COMPRESSED_VECTOR(S1, U1), COMPRESSED_VECTOR(S2, U2) >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      DEBUG_PRECONDITION_EQUAL(size(x), size(y));
      typename Iterate<RHS>::result_type r = iterate(y);
      while (r)
      {
	 add_element(x, r.index(), *r);
	 ++r;
      }
   }
};

// subtract

template <typename LHS, typename RHS>
struct SubtractExpression
{
};

template <typename LHS, typename RHS, typename S1, typename U1, typename S2, typename U2>
struct SubtractInterface<LHS, RHS, VECTOR_EXPRESSION(S1, U1), VECTOR_EXPRESSION(S2, U2) >
   : SubtractExpression<LHS, RHS>
{
};

template <typename LHS, typename RHS, typename S1, typename U1, typename S2, typename U2>
struct SubtractInterface<LHS, RHS, VECTOR_EXPRESSION(S1, U1), LOCAL_VECTOR(S2, U2) >
{
};

template <typename LHS, typename RHS, typename S1, typename U1, typename S2, typename U2>
struct SubtractInterface<LHS&, RHS, DENSE_VECTOR(S1, U1), DENSE_VECTOR(S2, U2) >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      DEBUG_PRECONDITION_EQUAL(size(x), size(y));
      iter_subtract(iterate(x), iterate(y));
   }
};

template <typename LHS, typename RHS, typename S1, typename U1, typename S2, typename U2>
struct SubtractInterface<LHS&, RHS, DENSE_VECTOR(S1, U1), COMPRESSED_VECTOR(S2, U2) >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      DEBUG_PRECONDITION_EQUAL(size(x), size(y));
      typename Iterate<RHS>::result_type r = iterate(y);
      typename Iterate<LHS&>::result_type l = iterate(x);
      while (r)
      {
	 l[r.index()] -= *r;
	 ++r;
      }
   }
};

template <typename LHS, typename RHS, typename S1, typename U1, typename S2, typename U2>
struct SubtractInterface<LHS&, RHS, COMPRESSED_VECTOR(S1, U1), COMPRESSED_VECTOR(S2, U2) >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      PRECONDITION_EQUAL(size(x), size(y));
      typename Iterate<RHS>::result_type r = iterate(y);
      while (r)
      {
	 subtract_element(x, r.index(), *r);
	 ++r;
      }
   }
};

// multiply

// TODO: we need to define this for all vectors, otherwise it would revert to the default
// version which delegates to operator*=, and that could cause recursion.
// But, we could patch this by making it empty if the nested multiply does not exist.
template <typename LHS, typename RHS, typename S1, typename U1, typename RHSi>
struct MultiplyInterface<LHS&, RHS, LOCAL_VECTOR(S1, U1), AnyScalar<RHSi> >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef RHS const& second_argument_type;
   void operator()(LHS& x, RHS const& y) const
   {
      typename Iterate<LHS&>::result_type i = iterate(x);
      while (i)
      {
         multiply(*i, y);
         ++i;
      }
   }
};

//
// defaults for all vectors - delegate to member functions
//

template <typename T, typename Enable = void>
struct TryIterate { };

template <typename T>
struct TryIterate<T&, typename boost::enable_if<
  boost::mpl::and_<exists<typename T::iterator>,
		   boost::mpl::not_<boost::is_const<T> > > >::type>
{
   typedef typename T::iterator result_type;
   typedef T& argument_type;
   result_type operator()(T& v) const { return v.iterate(); }
};

template <typename T>
struct TryIterate<T, 
		  typename boost::enable_if<exists<typename T::const_iterator> >::type>
{
   typedef typename T::const_iterator result_type;
   typedef T const& argument_type;
   result_type operator()(T const& v) const { return v.iterate(); }
};

template <typename T>
struct Iterate : TryIterate<T> { };

template <typename T, typename Tv, typename Ti>
struct SizeInterface<T, VECTOR_EXPRESSION(Tv, Ti)>
{
   typedef size_type result_type;
   typedef T const& argument_type;
   size_type operator()(T const& v) const { return v.size(); }
};

//
// defaults for dense vectors
//

template <typename T, typename Value, typename S, typename U>
struct SetElement<T&, Value, DENSE_VECTOR(S, U) >
{
   typedef void result_type;
   typedef T& first_argument_type;
   typedef size_type second_argument_type;
   typedef Value const& third_argument_type;
   result_type operator()(T& v, size_type n, Value const& x) const
   {
      v[n] = x;
   }
};

template <typename T, typename Value, typename S, typename U>
struct AddElement<T&, Value, DENSE_VECTOR(S, U) >
{
   typedef void result_type;
   typedef T& first_argument_type;
   typedef size_type second_argument_type;
   typedef Value const& third_argument_type;
   result_type operator()(T& v, size_type n, Value const& x) const
   {
      v[n] += x;
   }
};

template <typename T, typename Value, typename S, typename U>
struct SubtractElement<T&, Value, DENSE_VECTOR(S, U) >
{
   typedef void result_type;
   typedef T& first_argument_type;
   typedef size_type second_argument_type;
   typedef Value const& third_argument_type;
   result_type operator()(T& v, size_type n, Value const& x) const
   {
      v[n] -= x;
   }
};

template <typename T, typename S, typename U>
struct ZeroElement<T&, DENSE_VECTOR(S, U) >
{
   typedef void result_type;
   typedef T& first_argument_type;
   typedef size_type second_argument_type;
   result_type operator()(T& v, size_type n) const
   {
      v[n] = zero_or_die<S>();
   }
};

template <typename T, typename S, typename U>
struct ZeroAllInterface<T&, DENSE_VECTOR(S, U) >
{
   typedef void result_type;
   typedef T& argument_type;
   void operator()(T& v) const
   {
      iter_zero(iterate(v));
   }
};

template <typename T, typename S, typename U>
struct NNZInterface<T, DENSE_VECTOR(S, U) > : Size<T> { };

//
// defaults for sparse vectors, delegate to member functions
//

template <typename T, typename S, typename U>
struct NNZInterface<T, COMPRESSED_VECTOR(S, U)>
{
   typedef size_type result_type;
   typedef T argument_type;
   size_type operator()(T const& v) const { return v.nnz(); }
};

template <typename T, typename Value, typename S, typename U>
struct SetElement<T&, Value, INJECTIVE_VECTOR(S, U) >
{
   typedef void result_type;
   typedef T& first_argument_type;
   typedef size_type second_argument_type;
   typedef Value const& third_argument_type;
   result_type operator()(T& v, size_type n, Value const& x) const
   {
      v.set_element(n, x);
   }
};

template <typename T, typename Value, typename S, typename U>
struct AddElement<T&, Value, COMPRESSED_VECTOR(S, U) >
{
   typedef void result_type;
   typedef T& first_argument_type;
   typedef size_type second_argument_type;
   typedef Value const& third_argument_type;
   result_type operator()(T& v, size_type n, Value const& x) const
   {
      v.add_element(n, x);
   }
};

template <typename T, typename Value, typename S, typename U>
struct SubtractElement<T&, Value, COMPRESSED_VECTOR(S, U) >
{
   typedef void result_type;
   typedef T& first_argument_type;
   typedef size_type second_argument_type;
   typedef Value const& third_argument_type;
   result_type operator()(T& v, size_type n, Value const& x) const
   {
      v.subtract_element(n, x);
   }
};

template <typename T, typename S, typename U>
struct ZeroElement<T&, INJECTIVE_VECTOR(S, U) >
{
   typedef void result_type;
   typedef T& first_argument_type;
   typedef size_type second_argument_type;
   result_type operator()(T& v, size_type n) const
   {
      v.zero_element(n);
   }
};

template <typename T, typename S, typename U>
struct ZeroAllInterface<T&, VECTOR_EXPRESSION(S, U) >
{
   typedef void result_type;
   typedef T& argument_type;
   void operator()(T& v) const
   {
      v.zero_all();
   }
};

} // namespace LinearAlgebra

//#include "vectoroperations.cc"

#endif
