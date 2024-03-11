// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/matrixoperationsbase.h
//
// Copyright (C) 2004-2017 Ian McCulloch <ian@qusim.net>
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

#if !defined(MATRIXOPERATIONSBASE_H_YEWR4375678YOYOW)
#define MATRIXOPERATIONSBASE_H_YEWR4375678YOYOW

#include "matrixinterface.h"
#include "vectoroperationsbase.h"
#include "iteroperations2.h"
#include <boost/mpl/print.hpp>

namespace LinearAlgebra
{

// EvalExpression

template <typename T, typename V, typename U>
struct EvalExpression<T, Concepts::LocalMatrix<V, U>> : Identity<T> {};

template <typename T, typename V, typename U>
struct EvalExpression<T, Concepts::MatrixExpression<V, U>>
{
   typedef typename make_value<T>::type result_type;
   typedef T argument_type;
   result_type operator()(T const& x) const { return result_type(x); }
};

// Size1

template <typename T, typename TInterface = typename interface<T>::type>
struct Size1Interface { };

template <typename T, typename Enable = void>
struct Size1 : Size1Interface<T> {};

template <typename T>
inline
typename Size1<T>::result_type
size1(T const& x)
{
   return Size1<T>()(x);
}

// Size2

template <typename T, typename TInterface = typename interface<T>::type>
struct Size2Interface { };

template <typename T, typename Enable = void>
struct Size2 : Size2Interface<T> {};

template <typename T>
inline
typename Size2<T>::result_type
size2(T const& x)
{
   return Size2<T>()(x);
}

// Stride1 (Stride-matrix interface)

template <typename T, typename TInterface = typename interface<T>::type,
          typename Enable = void>
struct Stride1 { };

template <typename T>
inline
typename Stride1<T>::result_type
stride1(T const& x)
{
   return Stride1<T>()(x);
}

// Stride2 (Stride-matrix interface)

template <typename T, typename TInterface = typename interface<T>::type,
          typename Enable = void>
struct Stride2 { };

template <typename T>
inline
typename Stride2<T>::result_type
stride2(T const& x)
{
   return Stride2<T>()(x);
}

// Resize

template <typename T>
inline
typename boost::enable_if<is_matrix<T>, Resize<T&> >::type::result_type
resize(T& v, size_type r, size_type c)
{
   return Resize<T&>()(v, r, c);
}

template <typename T>
inline
typename boost::enable_if<is_matrix<T>, Resize<T&> >::type::result_type
try_resize(T& v, size_type r, size_type c)
{
   return Resize<T&>()(v, r, c);
}

template <typename T>
inline
typename boost::disable_if<is_defined<Resize<T&> >, void>::type
try_resize(T& v, size_type r, size_type c)
{
   PRECONDITION_EQUAL(size1(v), r);
   PRECONDITION_EQUAL(size2(v), c);
}

// Fill default implementation

template <typename T, typename V, typename Tv, typename Orient, typename Ti>
struct Fill<T&, V, Concepts::DenseMatrix<Tv, Orient, Ti>>
{
   typedef void result_type;
   typedef T& first_argument_type;
   typedef V second_argument_type;
   void operator()(T& v, V const& x) const
   {
      iter_fill_2(iterate(v), x);
   }
};


// SwapSortOrder

template <typename T, typename TInterface = typename interface<T>::type>
struct SwapSortOrderInterface {};

template <typename T>
struct SwapSortOrder : SwapSortOrderInterface<T> {};

template <typename T>
struct SwapSortOrderRef : SwapSortOrder<T&> {};

template <typename T>
inline
typename SwapSortOrder<T>::result_type
swap_sort_order(T const& x)
{
   return SwapSortOrder<T>()(x);
}

template <typename T>
inline
typename SwapSortOrder<T&>::result_type
swap_sort_order(T& x)
{
   return SwapSortOrder<T&>()(x);
}

//
// Unary operations
//

// most unary operations should follow the same approach as for vectors:

template <typename T, typename S, typename U>
struct NegateInterface<T, Concepts::MatrixExpression<S,U>> : TransformIfDefined<T, Negate<S> > {};

template <typename T, typename S, typename U>
struct RealInterface<T, Concepts::MatrixExpression<S,U>> : TransformIfDefined<T, Real<S> > {};

template <typename T, typename S, typename U>
struct RealInterface<T&, Concepts::MatrixExpression<S,U>> : TransformIfDefined<T&, Real<S&> > {};

template <typename T, typename S, typename U>
struct ImagInterface<T, Concepts::MatrixExpression<S,U>> : TransformIfDefined<T, Imag<S> > {};

template <typename T, typename S, typename U>
struct ImagInterface<T&, Concepts::MatrixExpression<S,U>> : TransformIfDefined<T&, Imag<S&> > {};

template <typename T, typename S, typename U>
struct ConjInterface<T, Concepts::MatrixExpression<S,U>> : TransformIfDefined<T, Conj<S> > {};

// Herm is conj(transpose)
template <typename T, typename S, typename U>
struct HermInterface<T, Concepts::MatrixExpression<S,U>>
   : UnaryComposer<Conj<typename Transpose<T>::result_type>,
                   Transpose<T> > {};

// Abs for matrices probably should NOT be defined.  Need to think about this.
template <typename T, typename S, typename U>
struct AbsInterface<T, Concepts::MatrixExpression<S,U>> : TransformIfDefined<T, Abs<S> > {};

// trace

template <typename T, typename Nested, typename Ti = typename interface<T>::type>
struct MatrixTrace {};

template <typename T, typename Tv, typename Ti>
struct TraceInterface<T, Concepts::AnyMatrix<Tv, Ti>> : MatrixTrace<T, Trace<Tv> > {};

// default implementation

// TODO: better implementation is possible in most cases
template <typename T, typename Nested, typename Tv, typename Ti>
struct MatrixTrace<T, Nested, Concepts::LocalMatrix<Tv, Ti>>
{
   typedef typename make_value<typename Nested::result_type>::type result_type;
   typedef T const& argument_type;
   typedef T const& first_argument_type;
   typedef Nested second_argument_type;

   result_type operator()(first_argument_type x, second_argument_type f) const
   {
      typedef typename make_value_with_zero<result_type>::type zval_type;
      zval_type Result = zero<zval_type>();

      typename const_iterator<T>::type I = iterate(x);
      while (I)
      {
         typename const_inner_iterator<T>::type J = iterate(I);
         while (J)
         {
            if (J.index1() == J.index2())
               add(Result, f(*J));

            ++J;
         }
         ++I;
      }
      return Result;
   }

   result_type operator()(argument_type x) const
   {
      return operator()(x, Nested());
   }
};

// NormFrob

#if 0
template <typename T, typename Sv, typename Si>
struct NormFrobSq<T, Concepts::MatrixExpression<Sv, Si>>
 : NormFrobSq<typename EvalExpression<T>::result_type> {};
#else
template <typename T, typename Sv, typename Si>
struct NormFrobSq<T, Concepts::MatrixExpression<Sv, Si>>
{
   typedef NormFrobSq<typename EvalExpression<T>::result_type> fwd;
   typedef typename fwd::result_type result_type;
   typedef T const& argument_type;
   result_type operator()(argument_type x) const
   {
      return fwd()(x);
   }
};
#endif

template <typename T, typename Sv, typename Si>
                                           struct NormFrobSq<T, Concepts::LocalMatrix<Sv, Si>>
{
   typedef typename make_value<typename NormFrobSq<Sv>::result_type>::type result_type;
   typedef T const& argument_type;
   result_type operator()(argument_type x) const
   {
      // find the first non-zero matrix element
      typename const_iterator<T>::type I = iterate(x);
      if (!I) return zero_or_die<result_type>();
      typename const_inner_iterator<T>::type J = iterate(I);
      while (!J)
      {
         ++I;
         if (!I) return zero_or_die<result_type>();
         J = iterate(I);
      }
      result_type Acc = norm_frob_sq(*J);
      ++J;
      while (J)
      {
         Acc += norm_frob_sq(*J);
         ++J;
      }
      ++I;
      while (I)
      {
         J = iterate(I);
         while (J)
         {
            Acc += norm_frob_sq(*J);
            ++J;
         }
         ++I;
      }
      return Acc;
   }
};

#if 0
// This is not the infinity norm!  It is also buggy, if the first row is empty...

template <typename T, typename Sv, typename Si>
struct NormInf<T, Concepts::LocalMatrix<Sv, Si>>
{
   typedef typename make_value<typename NormInf<Sv>::result_type>::type result_type;
   typedef T const& argument_type;
   result_type operator()(argument_type x) const
   {
      typename const_iterator<T>::type I = iterate(x);
      if (!I) return zero_or_die<result_type>();
      typename const_inner_iterator<T>::type J = iterate(I);
      if (!J) return zero_or_die<result_type>();
      result_type Max = norm_inf(*J);
      ++J;
      while (J)
      {
         if( Max < norm_inf(*J) )  Max = norm_inf(*J);
         ++J;
      }
      ++I;
      while (I)
      {
         J = iterate(I);
         while (J)
         {
            if( Max < norm_inf(*J) )  Max = norm_inf(*J);
            ++J;
         }
         ++I;
      }
      return Max;
   }
};
#endif

// inner_prod

template <typename S, typename T,
          typename Nested = InnerProd<typename interface<S>::value_type,
                                      typename interface<T>::value_type>,
          typename SInterface = typename interface<S>::type,
          typename TInterface = typename interface<T>::type,
          typename Enable = void>
struct MatrixInnerProd {};

// overload of inner_prod for 3-argument version
template <typename S, typename T, typename Func>
inline
typename MatrixInnerProd<S, T, Func>::result_type
inner_prod(S const& x, T const& y, Func const& f)
{
   return MatrixInnerProd<S, T, Func>()(x,y,f);
}

// forward InnerProd for matrix arguments to MatrixInnerProd
template <typename T, typename U, typename Tv, typename Ti, typename Uv, typename Ui>
struct InnerProdInterface<T, U, Concepts::AnyMatrix<Tv, Ti>,
                          Concepts::AnyMatrix<Uv, Ui>> : MatrixInnerProd<T, U> {};

template <typename S, typename T, typename Func,
          typename Sv, typename Si, typename Tv, typename Ti>
struct MatrixInnerProd<S, T, Func, Concepts::MatrixExpression<Sv, Si>, Concepts::MatrixExpression<Tv, Ti>>
   : MatrixInnerProd<typename EvalExpression<S>::result_type,
                     typename EvalExpression<T>::result_type,
                     Func> {};

template <typename S, typename T, typename Func,
          typename Sv, typename Sorient, typename Si,
          typename Tv, typename Torient, typename Ti>
struct MatrixInnerProd<S, T, Func, Concepts::DenseMatrix<Sv, Sorient, Si>,
                       Concepts::DenseMatrix<Tv, Torient, Ti>>
{
   typedef typename make_value_with_zero<typename Func::result_type>::type result_type;
   typedef S const& first_argument_type;
   typedef T const& second_argument_type;

   MatrixInnerProd() {}
   MatrixInnerProd(Func const& f) : f_(f) {}

   result_type operator()(S const& x, T const& y) const
   { return iter_inner_prod(iterate(x), iterate(swap_sort_order(y)),
                            VectorInnerProd<
                            typename Iterate<S>::result_type::reference,
                            typename Iterate<typename SwapSortOrder<T>::result_type>::result_type::reference,
                            Func>(f_)); }
   Func f_;
};

template <typename S, typename T, typename Func, typename Orient,
          typename Sv, typename Si, typename Tv, typename Ti>
struct MatrixInnerProd<S, T, Func, Concepts::DenseMatrix<Sv, Orient, Si>,
                       Concepts::DenseMatrix<Tv, Orient, Ti>>
{
   typedef typename make_value_with_zero<typename Func::result_type>::type result_type;
   typedef S const& first_argument_type;
   typedef T const& second_argument_type;

   MatrixInnerProd() {}
   MatrixInnerProd(Func const& f) : f_(f) {}

   result_type operator()(S const& x, T const& y) const
   { return iter_inner_prod(iterate(x), iterate(y),
                            VectorInnerProd<
                            typename Iterate<S>::result_type::reference,
                            typename Iterate<T>::result_type::reference,
                            Func>(f_)); }

   Func f_;
};

// This function looks completely bogus.  When did it get added?
template <typename S, typename T, typename Func, typename Orient,
          typename Sv, typename Si, typename Tv, typename Ti>
struct MatrixInnerProd<S, T, Func, Concepts::ContiguousMatrix<Sv, Orient, Si>,
                       Concepts::ContiguousMatrix<Tv, Orient, Ti>>
{
   typedef typename make_value_with_zero<typename Func::result_type>::type result_type;
   typedef S const& first_argument_type;
   typedef T const& second_argument_type;

   MatrixInnerProd() {}
   MatrixInnerProd(Func const& f) : f_(f) {}

   result_type operator()(S const& x, T const& y) const
   {
      return iter_inner_prod(VecPtrIterator<Sv const>(data(x), size1(x)*size2(x), 0),
                             VecPtrIterator<Tv const>(data(y), size1(y)*size2(y), 0),
                             f_);
   }
   Func f_;
};

template <typename S, typename T, typename Func, typename Orient,
          typename Sv, typename Si, typename Tv, typename Ti>
struct MatrixInnerProd<S, T, Func,
                       Concepts::CompressedOuterMatrix<Sv, Orient, Si>,
                       Concepts::CompressedOuterMatrix<Tv, Orient, Ti>>
{
   typedef typename make_value_with_zero<typename Func::result_type>::type result_type;
   typedef S const& first_argument_type;
   typedef T const& second_argument_type;

   MatrixInnerProd() {}
   MatrixInnerProd(Func const& f) : f_(f) {}

   result_type operator()(S const& x, T const& y) const
   { return iter_inner_prod(iterate(x), iterate(y),
                            VectorInnerProd<
                            typename Iterate<S>::result_type::reference,
                            typename Iterate<T>::result_type::reference,
                            Func>(f_)); }

   Func f_;
};

// parallel_prod
//
// for a matrix, this is equivalent to trace(A * transpose(B))

template <typename S, typename T, typename Nested,
          typename Sv, typename Si,
          typename Tv, typename Ti>
struct ParallelProdInterface<S, T, Nested, Concepts::AnyMatrix<Sv, Si>, Concepts::AnyMatrix<Tv, Ti>>
   : MatrixInnerProd<S, T, Nested> { };

// equal

template <typename T, typename U, typename TolType,
          typename Tv, typename Ti,
          typename Uv, typename Ui>
struct EqualInterface<T, U, TolType, Concepts::AnyMatrix<Tv,Ti>, Concepts::AnyMatrix<Uv, Ui>>
{
   typedef bool result_type;
   typedef T const& first_argument_type;
   typedef U const& second_argument_type;

   EqualInterface(TolType Tol = default_tolerance()) : Tol_(Tol) {}

   bool operator()(T const& x, U const& y) const
   {
      return (size1(x) == size1(y)) && (size2(x) == size2(y)) &&
         norm_frob(x-y) <= Tol_ + 2 * std::numeric_limits<TolType>::epsilon() *
         (norm_frob(x) + norm_frob(y));
   }

   private:
      TolType Tol_;
};

// equal_to

// TODO: equal_to for sparse matrices

template <typename S, typename T, typename Orient1,
          typename V1, typename U1,
          typename Orient2, typename V2, typename U2>
struct EqualToInterface<S, T, Concepts::DenseMatrix<V1, Orient1, U1>,
                        Concepts::DenseMatrix<V2, Orient2, U2>>
{
   typedef bool result_type;
   typedef S first_argument_type;
   typedef T second_argument_type;

   result_type operator()(S const& x, T const& y) const
   {
      return iter_equal_to(iterate(x), iterate(swap_sort_order(y)));
   }
};

template <typename S, typename T, typename Orient,
          typename V1, typename U1,
          typename V2, typename U2>
struct EqualToInterface<S, T, Concepts::DenseMatrix<V1, Orient, U1>,
                        Concepts::DenseMatrix<V2, Orient, U2>>
{
   typedef bool result_type;
   typedef S first_argument_type;
   typedef T second_argument_type;

   result_type operator()(S const& x, T const& y) const
   {
      return iter_equal_to(iterate(x), iterate(y));
   }
};

// multiplication

// matrix * scalar

template <typename S, typename T,
          typename NestedMultiply = Multiplication<typename interface<S>::value_type, T>,
          typename SInterface = typename interface<S>::type,
          typename Enable = void>
struct MatrixScalarMultiplication
{
};

template <typename S, typename T, typename Mult, typename Enable = void>
struct MatrixScalarMultiplicationImpl {};

template <typename S, typename T, typename Mult>
struct MatrixScalarMultiplicationImpl<S, T, Mult,
                                      typename boost::enable_if<is_defined<Mult> >::type>
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

template <typename S, typename T, typename SV, typename U, typename Mult>
struct MatrixScalarMultiplication<S, T, Mult, Concepts::AnyMatrix<SV, U>>
   : MatrixScalarMultiplicationImpl<S, T, Mult> {};

// we forward to MatrixScalarMultiplication<S, typename interface<T>::value_type>
// rather than using T directly, to handle the case where T is a ScalarProxy.

template <typename S, typename T, typename SV, typename SI>
struct MultiplicationInterface<S, T, Concepts::AnyMatrix<SV, SI>, AnyScalar<T> >
   : MatrixScalarMultiplication<S, typename interface<T>::value_type>
{
};

// scalar * matrix

template <typename S, typename T,
          typename NestedMultiply = Multiplication<S, typename interface<T>::value_type>,
          typename TInterface = typename interface<T>::type,
          typename Enable = void>
struct ScalarMatrixMultiplication
{
};

// backwards-compatability : better would be to implement scalar() function
template <typename S, typename T, typename Mult>
inline
typename ScalarMatrixMultiplication<S, T>::result_type
left_scalar_prod(S const& x, T const& y)
{
   return  ScalarMatrixMultiplication<S, T>()(x,y);
}

template <typename S, typename T, typename Mult>
inline
typename ScalarMatrixMultiplication<S, T, Mult>::result_type
left_scalar_prod(S const& x, T const& y, Mult const& f)
{
   return  ScalarMatrixMultiplication<S, T, Mult>()(x,y,f);
}

template <typename S, typename T, typename Mult, typename Enable = void>
struct ScalarMatrixMultiplicationImpl {};

template <typename S, typename T, typename Mult>
struct ScalarMatrixMultiplicationImpl<S, T, Mult,
                                      typename boost::enable_if<is_defined<Mult> >::type>
{
   typedef typename BindFirst<Mult>::result_type Func;
   typedef S first_argument_type;
   typedef T second_argument_type;
   typedef Mult third_argument_type;
   typedef typename Transform<T, Func>::result_type result_type;

   result_type operator()(S const& x, T const& v) const
   {
      return transform(v, bind_first(Mult(), x));
   }

   result_type operator()(S const& x, T const& v, Mult const& f) const
   {
      return transform(v, bind_first(f, x));
   }
};

template <typename S, typename T, typename TV, typename U, typename Mult>
struct ScalarMatrixMultiplication<S, T, Mult, Concepts::AnyMatrix<TV, U>>
   : ScalarMatrixMultiplicationImpl<S, T, Mult> {};

template <typename S, typename T, typename TV, typename TI>
struct MultiplicationInterface<S, T, AnyScalar<S>, Concepts::AnyMatrix<TV, TI>>
   : ScalarMatrixMultiplication<typename interface<S>::value_type, T> {};

// matrix-matrix multiply - declared here but defined in matrixmultiplication.h

template <typename S, typename T,
          typename NestedMultiply = Multiplication<typename interface<S>::value_type,
                                                   typename interface<T>::value_type>,
          typename Si = typename interface<S>::type,
          typename Ti = typename interface<T>::type>
struct MatrixMatrixMultiplication
{
};

template <typename S, typename T,
          typename Sv, typename Si,
          typename Tv, typename Ti>
struct MultiplicationInterface<S, T, Concepts::AnyMatrix<Sv,Si>, Concepts::AnyMatrix<Tv,Ti>>
   : MatrixMatrixMultiplication<S, T> {};

// transpose is defined in matrixtranspose.h

// get/set element interface for matrices

template <typename T, typename TInterface = typename interface<T>::type>
struct ZeroMatrixElementInterface {};

template <typename T, typename Enable = void>
struct ZeroMatrixElement : ZeroMatrixElementInterface<T> {};

template <typename T>
inline
typename ZeroMatrixElement<T&>::result_type
zero_element(T& m, size_type i, size_type j)
{
   return ZeroMatrixElement<T&>()(m, i, j);
}

// DirectSum

template <typename S, typename T,
          typename Nested = DirectSum<typename interface<S>::value_type,
                                      typename interface<T>::value_type>,
          typename SInterface = typename interface<S>::type,
          typename TInterface = typename interface<T>::type>
struct MatrixDirectSum {};

template <typename S, typename T, typename Sv, typename Si, typename Tv, typename Ti>
struct DirectSumInterface<S, T, Concepts::AnyMatrix<Sv, Si>, Concepts::AnyMatrix<Tv, Ti>>
   : MatrixDirectSum<S, T> {};

// DirectProduct

template <typename S, typename T,
          typename Nested = DirectProduct<typename interface<S>::value_type,
                                          typename interface<T>::value_type>,
          typename SInterface = typename interface<S>::type,
          typename TInterface = typename interface<T>::type>
struct MatrixDirectProduct {};

template <typename S, typename T, typename Sv, typename Si, typename Tv, typename Ti>
struct DirectProductInterface<S, T, Concepts::AnyMatrix<Sv, Si>, Concepts::AnyMatrix<Tv, Ti>>
   : MatrixDirectProduct<S, T> {};

template <typename T1, typename T2,
          typename T1i,
          typename T2v, typename T2i>
struct DirectProductInterface<T1, T2, AnyScalar<T1i>, Concepts::AnyMatrix<T2v, T2i>>
   : ScalarMatrixMultiplication<T1, T2> {};

template <typename T1, typename T2,
          typename T1v, typename T1i,
          typename T2i>
struct DirectProductInterface<T1, T2, Concepts::AnyMatrix<T1v, T1i>, AnyScalar<T2i>>
   : MatrixScalarMultiplication<T1, T2> {};

template <typename S, typename T, typename Nested>
inline
typename MatrixDirectProduct<S, T, Nested>::result_type
direct_product(S const& x, T const& y, Nested const& f)
{
   return MatrixDirectProduct<S, T, Nested>()(x,y, f);
}

// get_element

template <typename T, typename TInterface = typename interface<T>::type>
struct GetMatrixElementInterface {};

template <typename T, typename Enable = void>
struct GetMatrixElement : GetMatrixElementInterface<T> {};

template <typename T>
inline
typename boost::disable_if<is_mutable_proxy<T>, GetMatrixElement<T> >::type::result_type
get_element(T const& m, size_type i, size_type j)
{
   return GetMatrixElement<T>()(m, i, j);
}

template <typename T>
inline
typename GetMatrixElement<T&>::result_type
get_element(T& m, size_type i, size_type j)
{
   return GetMatrixElement<T&>()(m, i, j);
}

template <typename T>
inline
typename boost::enable_if<is_mutable_proxy<T>, GetMatrixElement<T&> >::type::result_type
get_element(T const& m, size_type i, size_type j)
{
   return GetMatrixElement<T&>()(const_cast<T&>(m), i, j);
}


// default implementation

template <typename T, typename Enable = void>
struct GetMatrixElementDefault {};

template <typename T>
struct GetMatrixElementDefault<T, typename boost::enable_if<exists<typename T::const_reference> >::type>
{
   typedef typename T::const_reference result_type;
   typedef T const& first_argument_type;
   typedef size_type second_argument_type;
   typedef size_type third_argument_type;
   result_type operator()(T const& m, size_type i, size_type j) const
   {
      return m(i,j);
   }
};

template <typename T>
struct GetMatrixElementDefault<T&, typename boost::enable_if<exists<typename T::reference> >::type>
{
   typedef typename T::reference result_type;
   typedef T& first_argument_type;
   typedef size_type second_argument_type;
   typedef size_type third_argument_type;
   result_type operator()(T& m, size_type i, size_type j) const
   {
      return m(i,j);
   }
};

template <typename T, typename Tv, typename Ti>
struct GetMatrixElementInterface<T, Concepts::AnyMatrix<Tv, Ti>> : GetMatrixElementDefault<T> {};

// set_element

template <typename T, typename V, typename TInterface = typename interface<T>::type>
struct SetMatrixElementInterface {};

template <typename T, typename V, typename Enable = void>
struct SetMatrixElement : SetMatrixElementInterface<T, V> {};

template <typename T, typename V>
inline
typename SetMatrixElement<T&, V>::result_type
set_element(T& m, size_type i, size_type j, V const& x)
{
   return SetMatrixElement<T&,V>()(m, i, j, x);
}

// set_element_lock

template <typename T, typename V, typename TInterface = typename interface<T>::type>
struct SetMatrixElementLockInterface {};

template <typename T, typename V, typename Enable = void>
struct SetMatrixElementLock : SetMatrixElementLockInterface<T, V> {};

template <typename T, typename V>
inline
typename SetMatrixElementLock<T&, V>::result_type
set_element_lock(T& m, size_type i, size_type j, V const& x)
{
   return SetMatrixElementLock<T&,V>()(m, i, j, x);
}

#if 0
// TODO: ZeroMatrixElementLock<> function
template <typename T, typename V>
inline
typename SetMatrixElementLock<T&, V>::result_type
set_element(T& m, size_type i, size_type j, value_with_zero<V> const& x)
{
   if (x.is_zero())
      ZeroMatrixElementLock<T&>(m,i,j);
   else
      SetMatrixElementLock<T&,V>()(m, i, j, x.get());
}
#endif

// TODO: this simply forwards to SetMatrixElement, it should instead
// forward to a new function SetNewMatrixElement
template <typename T, typename V>
inline
typename SetMatrixElement<T&, V>::result_type
set_new_element(T& m, size_type i, size_type j, V const& x)
{
   return SetMatrixElement<T&,V>()(m, i, j, x);
}

// TODO: this simply forwards to SetMatrixElement, it should instead
// forward to a new function SetNewMatrixElement
template <typename T, typename V>
inline
typename SetMatrixElement<T&, V>::result_type
set_new_element(T& m, size_type i, size_type j, value_with_zero<V> const& x)
{
   if (x.is_zero())
      ZeroMatrixElement<T&>(m,i,j);
   else
      SetMatrixElement<T&,V>()(m, i, j, x.get());
}

template <typename T, typename V>
inline
typename SetMatrixElement<T&, V>::result_type
set_element_check_if_zero(T& m, size_type i, size_type j, V const& x)
{
   if (is_zero(x))
      ZeroMatrixElement<T&>(m,i,j);
   else
      SetMatrixElement<T&,V>()(m, i, j, x);
}

template <typename T, typename V>
inline
typename SetMatrixElement<T&, V>::result_type
set_element_check_if_zero(T& m, size_type i, size_type j, value_with_zero<V> const& x)
{
   if (is_zero(x))
      ZeroMatrixElement<T&>(m,i,j);
   else
      SetMatrixElement<T&,V>()(m, i, j, x.get());
}

// default implementation

template <typename T, typename V, typename Enable = void>
struct SetMatrixElementDefault {};

template <typename T, typename V>
struct SetMatrixElementDefault<
   T&
 , V
 , typename boost::enable_if<is_defined<GetMatrixElement<T&> > >::type
>
{
   typedef void result_type;
   typedef T& first_argument_type;
   typedef size_type second_argument_type;
   typedef size_type third_argument_type;
   typedef V const& fourth_argument_type;
   void operator()(T& m, size_type i, size_type j, V const& x) const
   {
      get_element(m,i,j) = x;
   }
};

template <typename T, typename V, typename Tv, typename Ti>
struct SetMatrixElementInterface<T, V, Concepts::AnyMatrix<Tv, Ti>> : SetMatrixElementDefault<T, V> {};

// add_element

template <typename T, typename V, typename TInterface = typename interface<T>::type>
struct AddMatrixElementInterface {};

template <typename T, typename V, typename Enable = void>
struct AddMatrixElement : AddMatrixElementInterface<T, V> {};

template <typename T, typename V>
inline
typename AddMatrixElement<T&, V>::result_type
add_element(T& m, size_type i, size_type j, V const& x)
{
   return AddMatrixElement<T&,V>()(m, i, j, x);
}

template <typename T, typename V>
inline
typename AddMatrixElement<T&, V>::result_type
add_element(T& m, size_type i, size_type j, value_with_zero<V> const& x)
{
   if (!x.is_zero())  AddMatrixElement<T&,V>()(m, i, j, x);
}

template <typename T, typename V>
inline
typename AddMatrixElement<T&, V>::result_type
add_element_check_if_zero(T& m, size_type i, size_type j, V const& x)
{
   if (!is_zero(x)) AddMatrixElement<T&,V>()(m, i, j, x);
}

template <typename T, typename V>
inline
typename AddMatrixElement<T&, V>::result_type
add_element_check_if_zero(T& m, size_type i, size_type j, value_with_zero<V> const& x)
{
   if (!is_zero(x))  AddMatrixElement<T&,V>()(m, i, j, x);
}

// subtract_element

template <typename T, typename V, typename TInterface = typename interface<T>::type>
struct SubtractMatrixElementInterface {};

template <typename T, typename V, typename Enable = void>
struct SubtractMatrixElement : SubtractMatrixElementInterface<T, V> {};

template <typename T, typename V>
inline
typename SubtractMatrixElement<T&, V>::result_type
subtract_element(T& m, size_type i, size_type j, V const& x)
{
   return SubtractMatrixElement<T&,V>()(m, i, j, x);
}

template <typename T, typename V>
inline
typename SubtractMatrixElement<T&, V>::result_type
subtract_element(T& m, size_type i, size_type j, value_with_zero<V> const& x)
{
   if (!x.is_zero())  SubtractMatrixElement<T&,V>()(m, i, j, x);
}

template <typename T, typename V>
inline
typename SubtractMatrixElement<T&, V>::result_type
subtract_element_check_if_zero(T& m, size_type i, size_type j, V const& x)
{
   if (!is_zero(x)) SubtractMatrixElement<T&,V>()(m, i, j, x);
}

template <typename T, typename V>
inline
typename SubtractMatrixElement<T&, V>::result_type
subtract_element_check_if_zero(T& m, size_type i, size_type j, value_with_zero<V> const& x)
{
   if (!is_zero(x)) SubtractMatrixElement<T&,V>()(m, i, j, x);
}

// multiply_element

template <typename T, typename V, typename TInterface = typename interface<T>::type>
struct MultiplyMatrixElementInterface {};

template <typename T, typename V, typename Enable = void>
struct MultiplyMatrixElement : MultiplyMatrixElementInterface<T, V> {};

template <typename T, typename V>
inline
typename SubtractMatrixElement<T&, V>::result_type
multiply_element(T& m, size_type i, size_type j, V const& x)
{
   return MultiplyMatrixElement<T&,V>()(m, i, j, x);
}

template <typename T, typename V>
inline
typename MultiplyMatrixElement<T&, V>::result_type
multiply_element(T& m, size_type i, size_type j, value_with_zero<V> const& x)
{
   if (x.is_zero())
      ZeroMatrixElement<T>(m,i,j);
   else
      MultiplyMatrixElement<T&,V>()(m, i, j, x);
}

template <typename T, typename V>
inline
typename SubtractMatrixElement<T&, V>::result_type
multiply_element_check_if_zero(T& m, size_type i, size_type j, V const& x)
{
   if (is_zero(x))
      ZeroMatrixElement<T>(m,i,j);
   else
      MultiplyMatrixElement<T&,V>()(m, i, j, x);
}

template <typename T, typename V>
inline
typename MultiplyMatrixElement<T&, V>::result_type
multiply_element_check_if_zero(T& m, size_type i, size_type j, value_with_zero<V> const& x)
{
   if (is_zero(x))
      ZeroMatrixElement<T>(m,i,j);
   else
      MultiplyMatrixElement<T&,V>()(m, i, j, x);
}

// default implementation

template <typename T, typename V, typename Enable = void>
struct AddMatrixElementDefault {};

template <typename T, typename V>
struct AddMatrixElementDefault<
   T&
 , V
 , typename boost::enable_if<is_defined<GetMatrixElement<T&> > >::type
>
{
   typedef void result_type;
   typedef T& first_argument_type;
   typedef size_type second_argument_type;
   typedef size_type third_argument_type;
   typedef V const& fourth_argument_type;
   void operator()(T& m, size_type i, size_type j, V const& x) const
   {
      get_element(m,i,j) += x;
   }
};

template <typename T, typename V, typename Tv, typename Ti>
struct AddMatrixElementInterface<T, V, Concepts::AnyMatrix<Tv, Ti>> : AddMatrixElementDefault<T, V> {};

// default implementation

template <typename T, typename V, typename Enable = void>
struct SubtractMatrixElementDefault {};

template <typename T, typename V>
struct SubtractMatrixElementDefault<
   T&
 , V
 , typename boost::enable_if<is_defined<GetMatrixElement<T&> > >::type
>
{
   typedef void result_type;
   typedef T& first_argument_type;
   typedef size_type second_argument_type;
   typedef size_type third_argument_type;
   typedef V const& fourth_argument_type;
   void operator()(T& m, size_type i, size_type j, V const& x) const
   {
      get_element(m,i,j) -= x;
   }
};

template <typename T, typename V, typename Tv, typename Ti>
struct SubtractMatrixElementInterface<T, V, Concepts::AnyMatrix<Tv, Ti>>
   : SubtractMatrixElementDefault<T, V> {};

// default implementation

template <typename T, typename Tv = typename interface<T>::value_type, typename Enable = void>
struct ZeroMatrixElementDenseDefault {};

template <typename T, typename Tv>
struct ZeroMatrixElementDenseDefault<T, Tv, typename boost::enable_if<has_zero<Tv> >::type>
{
   typedef void result_type;
   typedef T& first_argument_type;
   typedef size_type second_argument_type;
   typedef size_type third_argument_type;
   result_type operator()(T& m, size_type i, size_type j)
   {
      set_element(m,i,j,zero<Tv>());
   }
};

template <typename T, typename Tv, typename TOrient, typename Ti>
struct ZeroMatrixElementInterface<T, Concepts::DenseMatrix<Tv,TOrient,Ti>>
   : ZeroMatrixElementDenseDefault<T, Tv> {};

// vector_view - views a matrix as a vector of vectors.
// only applicable to matrices that have a definite major ordering.

template <typename T, typename Ti = typename interface<T>::type>
struct VectorViewInterface {};

template <typename T, typename Enable = void>
struct VectorView : VectorViewInterface<T> {};

template <typename T>
inline
typename VectorView<T>::result_type
vector_view(T const& x)
{
   return VectorView<T>()(x);
}

template <typename T>
inline
typename VectorView<T&>::result_type
vector_view(T& x)
{
   return VectorView<T&>()(x);
}

// flatten_rows and flatten_cols
// represent a matrix as a vector

template <typename T, typename Ti = typename interface<T>::type>
struct FlattenRowsInterface {};

template <typename T>
struct FlattenRows : FlattenRowsInterface<T> {};

template <typename T>
inline
typename FlattenRows<T>::result_type
flatten_rows(T const& x)
{
   return FlattenRows<T>()(x);
}

template <typename T, typename Ti = typename interface<T>::type>
struct FlattenColsInterface {};

template <typename T>
struct FlattenCols : FlattenColsInterface<T> {};

template <typename T>
inline
typename FlattenCols<T>::result_type
flatten_cols(T const& x)
{
   return FlattenCols<T>()(x);
}

// default implementations

template <typename T, typename Tv, typename Ti>
struct FlattenRowsInterface<T, Concepts::ContiguousMatrix<Tv, RowMajor, Ti>>
{
   typedef VectorMemProxy<Tv const> result_type;
   result_type operator()(T const& m) const
   {
      return result_type(data(m), size1(m) * size2(m));
   }
};

template <typename T, typename Tv, typename Ti>
struct FlattenColsInterface<T, Concepts::ContiguousMatrix<Tv, ColMajor, Ti>>
{
   typedef VectorMemProxy<Tv const> result_type;
   result_type operator()(T const& m) const
   {
      return result_type(data(m), size1(m) * size2(m));
   }
};

// iterate_at

template <typename T, typename Ti = typename interface<T>::type>
struct MatrixIterateAtInterface {};

template <typename T, typename Enable = void>
struct MatrixIterateAt : MatrixIterateAtInterface<T> {};

template <typename T>
typename MatrixIterateAt<T>::result_type
iterate_at(T const& m, size_type i, size_type j)
{
   return MatrixIterateAt<T>()(m, i, j);
}

template <typename T>
typename MatrixIterateAt<T&>::result_type
iterate_at(T& m, size_type i, size_type j)
{
   return MatrixIterateAt<T&>()(m, i, j);
}

// default implementation for dense
template <typename T, typename Tv, typename Ti>
struct MatrixIterateAt<T, Concepts::DenseMatrix<Tv, RowMajor, Ti>>
{
   typedef typename const_inner_iterator<T>::type result_type;
   typedef T const& first_argument_type;
   typedef size_type second_argument_type;
   typedef size_type third_argument_type;

   result_type operator()(first_argument_type m,
                          second_argument_type i,
                          third_argument_type j) const
   {
      typename const_iterator<T>::type I = iterate(m);
      I += i;
      result_type J = iterate(I);
      J += j;
      return J;
   }
};

template <typename T, typename Tv, typename Ti>
struct MatrixIterateAt<T&, Concepts::DenseMatrix<Tv, RowMajor, Ti>>
{
   typedef typename inner_iterator<T>::type result_type;
   typedef T& first_argument_type;
   typedef size_type second_argument_type;
   typedef size_type third_argument_type;

   result_type operator()(first_argument_type m,
                          second_argument_type i,
                          third_argument_type j) const
   {
      typename iterator<T>::type I = iterate(m);
      I += i;
      result_type J = iterate(I);
      J += j;
      return J;
   }
};

template <typename T, typename Tv, typename Ti>
struct MatrixIterateAt<T, Concepts::DenseMatrix<Tv, ColMajor, Ti>>
{
   typedef typename const_inner_iterator<T>::type result_type;
   typedef T const& first_argument_type;
   typedef size_type second_argument_type;
   typedef size_type third_argument_type;

   result_type operator()(first_argument_type m,
                          second_argument_type i,
                          third_argument_type j) const
   {
      typename const_iterator<T>::type I = iterate(m);
      I += j;
      result_type J = iterate(I);
      J += i;
      return J;
   }
};

template <typename T, typename Tv, typename Ti>
struct MatrixIterateAt<T&, Concepts::DenseMatrix<Tv, ColMajor, Ti>>
{
   typedef typename inner_iterator<T>::type result_type;
   typedef T& first_argument_type;
   typedef size_type second_argument_type;
   typedef size_type third_argument_type;

   result_type operator()(first_argument_type m,
                          second_argument_type i,
                          third_argument_type j) const
   {
      typename iterator<T>::type I = iterate(m);
      I += j;
      result_type J = iterate(I);
      J += i;
      return J;
   }
};

// MatrixRow

template <typename T, typename TInterface = typename interface<T>::type>
struct MatrixRowInterface { };

template <typename T, typename Enable = void>
struct MatrixRow : MatrixRowInterface<T> { };

template <typename T>
inline
typename MatrixRow<T>::result_type
matrix_row(T const& x, size_type n)
{
   return MatrixRow<T>()(x,n);
}

template <typename T>
inline
typename MatrixRow<T&>::result_type
matrix_row(T& x, size_type n)
{
   return MatrixRow<T&>()(x,n);
}

// MatrixCol

template <typename T, typename TInterface = typename interface<T>::type>
struct MatrixColInterface { };

template <typename T, typename Enable = void>
struct MatrixCol : MatrixColInterface<T> { };

template <typename T>
inline
typename MatrixCol<T>::result_type
matrix_col(T const& x, size_type n)
{
   return MatrixCol<T>()(x,n);
}

template <typename T>
inline
typename MatrixCol<T&>::result_type
matrix_col(T& x, size_type n)
{
   return MatrixCol<T&>()(x,n);
}

// Index1

template <typename T, typename U,
          typename Ti = typename interface<T>::type,
          typename Ui = typename interface<U>::type>
struct Index1Interface {};

template <typename T, typename U, typename Enable = void>
struct Index1 : Index1Interface<T,U> {};

// Index2

template <typename T, typename U,
          typename Ti = typename interface<T>::type,
          typename Ui = typename interface<U>::type>
struct Index2Interface {};

template <typename T, typename U, typename Enable = void>
struct Index2 : Index2Interface<T,U> {};

// project1, project2

template <typename T, typename U,
          typename Ti = typename interface<T>::type,
          typename Ui = typename interface<U>::type>
struct Project1Interface {};

template <typename T, typename U, typename Enable = void>
struct Project1 : Project1Interface<T,U> {};

template <typename T, typename U,
          typename Tv, typename Ti,
          typename Uv, typename Ui>
struct Project1Interface<T, U, Concepts::AnyMatrix<Tv, Ti>, ANY_VECTOR(Uv, Ui)>
   : Index1<T,U> {};

template <typename T, typename U>
struct Project1<T, U, typename boost::enable_if<boost::is_convertible<U, size_type> >::type>
   : MatrixRow<T> {};

template <typename T, typename U>
inline
typename Project1<T, U>::result_type
project1(T const& m, U const& x)
{
   return Project1<T,U>()(m,x);
}

template <typename T, typename U>
inline
typename Project1<T&, U>::result_type
project1(T& m, U const& x)
{
   return Project1<T&,U>()(m,x);
}

// project2

template <typename T, typename U,
          typename Ti = typename interface<T>::type,
          typename Ui = typename interface<U>::type>
struct Project2Interface {};

template <typename T, typename U, typename Enable = void>
struct Project2 : Project2Interface<T,U> {};

template <typename T, typename U,
          typename Tv, typename Ti,
          typename Uv, typename Ui>
struct Project2Interface<T, U, Concepts::AnyMatrix<Tv, Ti>, ANY_VECTOR(Uv, Ui)>
   : Index2<T,U> {};

template <typename T, typename U>
struct Project2<T, U, typename boost::enable_if<boost::is_convertible<U, size_type> >::type>
   : MatrixCol<T> {};

template <typename T, typename U>
inline
typename Project2<T, U>::result_type
project2(T const& m, U const& x)
{
   return Project2<T,U>()(m,x);
}

template <typename T, typename U>
inline
typename Project2<T&, U>::result_type
project2(T& m, U const& x)
{
   return Project2<T&,U>()(m,x);
}

// ProjectMatrix

struct all_helper {};

typedef int all_helper::*all_t;

all_t const all = ((all_t)0);

//struct all_t {};

//extern all_t all;

template <typename T, typename U, typename V,
          typename Ti = typename interface<T>::type,
          typename Ui = typename interface<U>::type,
          typename Vi = typename interface<V>::type>
struct ProjectMatrixInterface {};

template <typename T, typename U, typename V, typename Enable = void>
struct ProjectMatrix : ProjectMatrixInterface<T,U,V> {};

template <typename T, typename U, typename V>
inline
typename ProjectMatrix<T,U,V>::result_type
project(T const& m, U const& u, V const& v)
{
   return ProjectMatrix<T,U,V>()(m,u,v);
}

template <typename T, typename U, typename V>
inline
typename ProjectMatrix<T&,U,V>::result_type
project(T& m, U const& u, V const& v)
{
   return ProjectMatrix<T&,U,V>()(m,u,v);
}

// forwards for 'all' sectioning

template <typename T, typename U>
struct ProjectMatrix<T, U, all_t,
                     typename boost::enable_if<is_defined<Project1<T,U> > >::type>
{
   typedef typename Project1<T,U>::result_type result_type;
   typedef T const& first_argument_type;
   typedef U const& second_argument_type;
   typedef all_t third_argument_type;
   result_type operator()(T const& m, U const& v, all_t) const
   {
      return Project1<T,U>()(m,v);
   }
};

template <typename T, typename U>
struct ProjectMatrix<T&, U, all_t,
                     typename boost::enable_if<is_defined<Project1<T&,U> > >::type>
{
   typedef typename Project1<T&,U>::result_type result_type;
   typedef T& first_argument_type;
   typedef U const& second_argument_type;
   typedef all_t third_argument_type;
   result_type operator()(T& m, U const& v, all_t) const
   {
      return Project1<T&,U>()(m,v);
   }
};

template <typename T, typename V>
struct ProjectMatrix<T, all_t, V,
                     typename boost::enable_if<is_defined<Project2<T,V> > >::type>
{
   typedef typename Project2<T,V>::result_type result_type;
   typedef T const& first_argument_type;
   typedef all_t second_argument_type;
   typedef V const& third_argument_type;
   result_type operator()(T const& m, all_t, V const& v) const
   {
      return Project2<T,V>()(m,v);
   }
};

template <typename T, typename V>
struct ProjectMatrix<T&, all_t, V,
                     typename boost::enable_if<is_defined<Project2<T&,V> > >::type>
{
   typedef typename Project2<T&,V>::result_type result_type;
   typedef T& first_argument_type;
   typedef all_t second_argument_type;
   typedef V const& third_argument_type;
   result_type operator()(T& m, all_t, V const& v) const
   {
      return Project2<T&,V>()(m,v);
   }
};

template <typename T, typename S1, typename S2>
struct ProjectMatrix<
   T
 , S1
 , S2
 , typename boost::enable_if<
      boost::mpl::and_<
         exists<typename GetMatrixElement<T>::result_type>
       , boost::mpl::not_<is_builtin<GetMatrixElement<T> > >
       , boost::is_arithmetic<S1>
       , boost::is_arithmetic<S2>
      >
   >::type
> : GetMatrixElement<T> {};

// defaults

template <typename T, typename U, typename Tv, typename Ti, typename Uv, typename Ui>
struct Index1Interface<T, U, Concepts::AnyMatrix<Tv, Ti>, ANY_VECTOR(Uv, Ui)>
{
   typedef T const& first_argument_type;
   typedef U const& second_argument_type;
   typedef typename ProjectMatrix<T, U, Range>::result_type result_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      return project(x, y, range(0,size2(x)));
   }
};

template <typename T, typename U, typename Tv, typename Ti, typename Uv, typename Ui>
struct Index1Interface<T&, U, Concepts::AnyMatrix<Tv, Ti>, ANY_VECTOR(Uv, Ui)>
{
   typedef T& first_argument_type;
   typedef U const& second_argument_type;
   typedef typename ProjectMatrix<T&, U, Range>::result_type result_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      return project(x, y, range(0,size2(x)));
   }
};

template <typename T, typename U, typename Tv, typename Ti, typename Uv, typename Ui>
struct Index2Interface<T, U, Concepts::AnyMatrix<Tv, Ti>, ANY_VECTOR(Uv, Ui)>
{
   typedef T const& first_argument_type;
   typedef U const& second_argument_type;
   typedef typename ProjectMatrix<T, Range, U>::result_type result_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      return project(x, range(0,size1(x)), y);
   }
};

template <typename T, typename U, typename Tv, typename Ti, typename Uv, typename Ui>
struct Index2Interface<T&, U, Concepts::AnyMatrix<Tv, Ti>, ANY_VECTOR(Uv, Ui)>
{
   typedef T& first_argument_type;
   typedef U const& second_argument_type;
   typedef typename ProjectMatrix<T&, Range, U>::result_type result_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      return project(x, range(0,size1(x)), y);
   }
};

//
// bracket operator
//

template <typename T, typename Arg1, typename Arg2, typename Enable = void>
struct MatrixBracket { };

template <typename T, typename Arg1, typename Arg2>
struct MatrixBracket<T, Arg1, Arg2> : ProjectMatrix<T, Arg1, Arg2> {};

//
// ZeroAll
//

template <typename T, typename S, typename Orient, typename U>
struct ZeroAllInterface<T&, Concepts::DenseMatrix<S, Orient, U>>
{
   typedef void result_type;
   typedef T& argument_type;
   void operator()(T& v) const
   {
      iter_zero(iterate(v));
   }
};

template <typename T, typename S, typename Orient, typename U>
struct ZeroAllInterface<T&, Concepts::CompressedOuterMatrix<S, Orient, U>>
{
   typedef void result_type;
   typedef T& argument_type;
   void operator()(T& v) const
   {
      iter_zero(iterate(v));
   }
};

template <typename T, typename S, typename U>
struct ZeroAllInterface<T&, Concepts::DiagonalMatrix<S, U>>
{
   typedef void result_type;
   typedef T& argument_type;
   void operator()(T& v) const
   {
      iter_zero(iterate(v.diagonal()));
   }
};

template <typename T, typename S, typename U>
struct ZeroAllInterface<T&, Concepts::ScalarMatrix<S, U>>
{
   typedef void result_type;
   typedef T& argument_type;
   void operator()(T& v) const
   {
      zero_all(v.value());
   }
};

// Max

template <typename T, typename Tv, typename Ti>
struct Max<T, Concepts::LocalMatrix<Tv, Ti>>
{
   typedef typename GetMatrixElement<T>::result_type result_type;
   typedef T const& argument_type;
   result_type operator()(argument_type v) const
   {
      return *iter_matrix_max(iterate(v));
   }
};

// Min

template <typename T, typename Tv, typename Ti>
struct Min<T, Concepts::LocalMatrix<Tv, Ti>>
{
   typedef typename GetMatrixElement<T>::result_type result_type;
   typedef T const& argument_type;
   result_type operator()(argument_type v) const
   {
      return *iter_matrix_min(iterate(v));
   }
};

// NormInf

template <typename T, typename Sv, typename Si>
struct NormInf<T, Concepts::MatrixExpression<Sv, Si>>
{
   typedef NormInf<typename EvalExpression<T>::result_type> fwd;
   typedef typename fwd::result_type result_type;
   typedef T const& argument_type;
   result_type operator()(argument_type x) const
   {
      return fwd()(x);
   }
};

//
// MatrixTranspose
//

template <typename M, typename NestedFunc, typename MInterface = typename interface<M>::type,
          typename Enable = void>
struct MatrixTransposeInterface {};

template <typename M, typename NestedFunc = Transpose<typename interface<M>::value_type>,
          typename Enable = void>
struct MatrixTranspose : MatrixTransposeInterface<M, NestedFunc> {};

template <typename M, typename Mv, typename Mi>
struct TransposeInterface<M, Concepts::AnyMatrix<Mv, Mi>> : MatrixTranspose<M> {};

template <typename M, typename Mv, typename Mi>
struct TransposeInterface<M&, Concepts::AnyMatrix<Mv, Mi>> : MatrixTranspose<M&> {};

template <typename M, typename NestedFunc>
inline
typename MatrixTranspose<M, NestedFunc>::result_type
transpose(M const& m, NestedFunc n)
{
   return MatrixTranspose<M, NestedFunc>()(m,n);
}

// specialization for diagonal_matrix, where transpose is the identity operation
// (assuming that the nested operation is trivial)

template <typename M, typename F, typename T, typename Ti>
struct MatrixTransposeInterface<M, F, Concepts::DiagonalMatrix<T, Ti>,
                                typename boost::enable_if<is_identity<F>>::type>
   : Identity<M> {};

// do we need a version with Identity<T&> ?  Probably not, do we really want to use trans() as an l-value?

//
// defaults
//

template <typename T, typename Tv, typename Ti>
struct Size1Interface<T, Concepts::AnyMatrix<Tv, Ti>>
{
   typedef size_type result_type;
   typedef T argument_type;
   size_type operator()(T const& m) const { return m.size1(); }
};

template <typename T, typename Tv, typename Ti>
struct Size2Interface<T, Concepts::AnyMatrix<Tv, Ti>>
{
   typedef size_type result_type;
   typedef T argument_type;
   size_type operator()(T const& m) const { return m.size2(); }
};

// stride

template <typename T, typename Ti>
struct Stride1<T, Ti>
{
   typedef size_type result_type;
   typedef T argument_type;
   size_type operator()(T const& m) const { return m.stride1(); }
};

template <typename T, typename Ti>
struct Stride2<T, Ti>
{
   typedef size_type result_type;
   typedef T argument_type;
   size_type operator()(T const& m) const { return m.stride2(); }
};

// data

template <typename T, typename Tv, typename To, typename Ti>
struct DataInterface<T, Concepts::StrideMatrix<Tv, To, Ti>>
{
   typedef Tv const* result_type;
   typedef T argument_type;
   result_type operator()(T const& x) const { return x.data(); }
};

template <typename T, typename Tv, typename To, typename Ti>
struct DataInterface<T&, Concepts::StrideMatrix<Tv, To, Ti>>
{
   typedef Tv* result_type;
   typedef T& argument_type;
   result_type operator()(T& x) const { return x.data(); }
};

// addition

template <typename S, typename T, typename Sv, typename Si, typename Tv, typename Ti>
struct AdditionInterface<S, T, Concepts::MatrixExpression<Sv, Si>, Concepts::MatrixExpression<Tv, Ti>>
   : BinaryTransform<S, T, Addition<Sv, Tv> > { };

// subtraction

template <typename S, typename T, typename Sv, typename Si, typename Tv, typename Ti>
struct SubtractionInterface<S, T, Concepts::MatrixExpression<Sv, Si>, Concepts::MatrixExpression<Tv, Ti>>
   : BinaryTransform<S, T, Subtraction<Sv, Tv> > { };

// Multiply

template <typename S, typename T, typename Sv, typename Si, typename Tv, typename Ti>
struct MultiplyInterface<S&, T, Concepts::MatrixExpression<Sv, Si>, Concepts::MatrixExpression<Tv, Ti>>
{
   typedef S& result_type;
   typedef S& first_argument_type;
   typedef T const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      assign_copy(x, x*y);
      return x;
   }
};

// element_prod

template <typename S, typename T>
struct ElementProd : BinaryTransform<S, T, Multiplication<typename interface<S>::value_type,
                                                          typename interface<T>::value_type> > {};

template <typename S, typename T>
typename ElementProd<S, T>::result_type
element_prod(S const& s, T const& t)
{
   return ElementProd<S, T>()(s,t);
}

//
// Assignment
//

template <typename LHS, typename RHS, typename S1, typename U1, typename S2, typename U2>
struct AssignInterface<LHS, RHS, Concepts::MatrixExpression<S1, U1>, Concepts::MatrixExpression<S2, U2>>
   : AssignExpression<LHS, RHS>
{
};

// generic/coordinate

template <typename LHS, typename RHS,
          typename S1, typename U1,
          typename S2, typename U2>
struct AssignInterface<LHS&, RHS, Concepts::LocalMatrix<S1, U1>, Concepts::LocalMatrix<S2, U2>>
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      try_resize(x, size1(y), size2(y));
      typedef typename const_iterator<RHS>::type outer_iterator;
      typedef typename const_iterator<outer_iterator>::type inner_iterator;
      zero_all(x);
      outer_iterator r = iterate(y);
      while (r)
      {
         inner_iterator s = iterate(r);
         while (s)
         {
            add_element(x, s.index1(), s.index2(), *s);
            ++s;
         }
         ++r;
      }
   }
};

template <typename LHS, typename RHS,
          typename S1, typename U1,
          typename S2, typename U2>
struct AssignInterface<LHS&, RHS, Concepts::LocalMatrix<S1, U1>, Concepts::InjectiveMatrix<S2, U2>>
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      try_resize(x, size1(y), size2(y));
      typedef typename const_iterator<RHS>::type outer_iterator;
      typedef typename const_iterator<outer_iterator>::type inner_iterator;
      zero_all(x);
      outer_iterator r = iterate(y);
      while (r)
      {
         inner_iterator s = iterate(r);
         while (s)
         {
            set_new_element(x, s.index1(), s.index2(), *s);
            ++s;
         }
         ++r;
      }
   }
};

// Diagonal matrices are a tricky one, if the value_type doesn't have a zero element
// but is zero-able (eg a DiagonalMatrix<Matrix<T>>,
// and the right hand side is sparse.  The implementation above doesn't quite do the right thing,
// because zero_all will zero all elements of the diagonal but leave them the same size.
// But it isn't clear what one WOULD expect in this case - if the right hand side contains every element
// then the result will be correct anyway, and if the right hand side has missing elements then it
// isn't clear what the diagonal should contain anyway.  It is basically a runtime error I think.

// base case covers mixed row/column major

template <typename LHS, typename RHS,
          typename S1, typename Orient1, typename U1,
          typename S2, typename Orient2, typename U2>
struct AssignInterface<LHS&, RHS,
                       Concepts::DenseMatrix<S1, Orient1, U1>,
                       Concepts::DenseMatrix<S2, Orient2, U2>>
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      return assign(x, swap_sort_order(y));
   }
};

template <typename LHS, typename RHS, typename S1, typename Orient1, typename U1,
          typename S2, typename Orient2, typename U2>
struct AssignInterface<LHS&, RHS,
                       Concepts::DenseMatrix<S1, Orient1, U1>,
                       Concepts::CompressedOuterMatrix<S2, Orient2, U2>>
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      return assign(swap_sort_order(x), y);
   }
};

// same-major matrices

template <typename LHS, typename RHS, typename Orient,
          typename S1, typename U1, typename S2, typename U2>
struct AssignInterface<LHS&, RHS,
                       Concepts::DenseMatrix<S1, Orient, U1>,
                       Concepts::DenseMatrix<S2, Orient, U2>>
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      try_resize(x, size1(y), size2(y));
      iter_assign_2(iterate(x), const_iterate(y));
   }
};

template <typename LHS, typename RHS, typename Orient,
          typename S1, typename U1, typename S2, typename U2>
struct AssignInterface<LHS&, RHS,
                       Concepts::CompressedOuterMatrix<S1, Orient, U1>,
                       Concepts::CompressedOuterMatrix<S2, Orient, U2>>
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      try_resize(x, size1(y), size2(y));
      zero_all(x);
      typename Iterate<RHS>::result_type r = iterate(y);
      while (r)
      {
         set_new_element(vector_view(x), r.index(), *r);
         ++r;
      }
   }
};

template <typename LHS, typename RHS, typename Orient,
          typename S1, typename U1, typename S2, typename U2>
struct AssignInterface<LHS&, RHS,
                       Concepts::DenseMatrix<S1, Orient, U1>,
                       Concepts::CompressedOuterMatrix<S2, Orient, U2>>
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      try_resize(x, size1(y), size2(y));
      typedef typename const_iterator<RHS>::type outer_iterator;
      typedef typename const_iterator<outer_iterator>::type inner_iterator;
      zero_all(x);
      outer_iterator r = iterate(y);
      while (r)
      {
         inner_iterator s = iterate(r);
         while (s)
         {
            add_element(x, s.index1(), s.index2(), *s);
            ++s;
         }
         ++r;
      }
   }
};

// addition

template <typename LHS, typename RHS, typename S1, typename U1, typename S2, typename U2>
struct AddInterface<LHS&, RHS, Concepts::MatrixExpression<S1, U1>, Concepts::MatrixExpression<S2, U2>>
   : AddExpression<LHS&, RHS>
{
};

template <typename LHS, typename RHS,
          typename S1, typename U1,
          typename S2, typename U2>
struct AddInterface<LHS&, RHS, Concepts::LocalMatrix<S1, U1>, Concepts::LocalMatrix<S2, U2>>
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      if (is_zero(x))
      {
	 assign(x, y);
	 return;
      }
      typedef typename const_iterator<RHS>::type outer_iterator;
      typedef typename const_iterator<outer_iterator>::type inner_iterator;
      outer_iterator r = iterate(y);
      while (r)
      {
         inner_iterator s = iterate(r);
         while (s)
         {
            add_element(x, s.index1(), s.index2(), *s);
            ++s;
         }
         ++r;
      }
   }
};

// base case covers mixed row/column major

template <typename LHS, typename RHS,
          typename S1, typename Orient1, typename U1,
          typename S2, typename Orient2, typename U2>
struct AddInterface<LHS&, RHS,
                    Concepts::DenseMatrix<S1, Orient1, U1>,
                    Concepts::DenseMatrix<S2, Orient2, U2>>
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      add(x, swap_sort_order(y));
   }
};

template <typename LHS, typename RHS, typename S1, typename Orient1, typename U1,
          typename S2, typename Orient2, typename U2>
struct AddInterface<LHS&, RHS,
                    Concepts::DenseMatrix<S1, Orient1, U1>,
                    Concepts::CompressedOuterMatrix<S2, Orient2, U2>>
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      add(swap_sort_order(x), y);
   }
};

// same-major

template <typename LHS, typename RHS, typename Orient,
          typename S1, typename U1, typename S2, typename U2>
struct AddInterface<LHS&, RHS,
                    Concepts::DenseMatrix<S1, Orient, U1>,
                    Concepts::DenseMatrix<S2, Orient, U2>>
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      PRECONDITION_EQUAL(size1(x), size1(y));
      PRECONDITION_EQUAL(size2(x), size2(y));
      iter_add(iterate(x), const_iterate(y));
   }
};

template <typename LHS, typename RHS, typename Orient,
          typename S1, typename U1, typename S2, typename U2>
struct AddInterface<LHS&, RHS,
                    Concepts::CompressedOuterMatrix<S1, Orient, U1>,
                    Concepts::CompressedOuterMatrix<S2, Orient, U2>>
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      PRECONDITION_EQUAL(size1(x), size1(y));
      PRECONDITION_EQUAL(size2(x), size2(y));
      typename Iterate<RHS>::result_type r = iterate(y);
      while (r)
      {
         add_element(vector_view(x), r.index(), *r);
         ++r;
      }
   }
};

template <typename LHS, typename RHS, typename Orient,
          typename S1, typename U1, typename S2, typename U2>
struct AddInterface<LHS&, RHS,
                    Concepts::DenseMatrix<S1, Orient, U1>,
                    Concepts::CompressedOuterMatrix<S2, Orient, U2>>
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      PRECONDITION_EQUAL(size1(x), size1(y));
      PRECONDITION_EQUAL(size2(x), size2(y));
      typedef typename const_iterator<RHS>::type outer_iterator;
      typedef typename const_iterator<outer_iterator>::type inner_iterator;
      outer_iterator r = iterate(y);
      while (r)
      {
         inner_iterator s = iterate(r);
         while (s)
         {
            add_element(x, s.index1(), s.index2(), *s);
            ++s;
         }
         ++r;
      }
   }
};

// subtract

template <typename LHS, typename RHS, typename S1, typename U1, typename S2, typename U2>
struct SubtractInterface<LHS, RHS, Concepts::MatrixExpression<S1, U1>, Concepts::MatrixExpression<S2, U2>>
   : SubtractExpression<LHS, RHS>
{
};

template <typename LHS, typename RHS,
          typename S1, typename U1,
          typename S2, typename U2>
struct SubtractInterface<LHS&, RHS, Concepts::LocalMatrix<S1, U1>, Concepts::LocalMatrix<S2, U2>>
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      typedef typename const_iterator<RHS>::type outer_iterator;
      typedef typename const_iterator<outer_iterator>::type inner_iterator;
      outer_iterator r = iterate(y);
      while (r)
      {
         inner_iterator s = iterate(r);
         while (s)
         {
            subtract_element(x, s.index1(), s.index2(), *s);
            ++s;
         }
         ++r;
      }
   }
};

// base case covers mixed row/column major

template <typename LHS, typename RHS,
          typename S1, typename Orient1, typename U1,
          typename S2, typename Orient2, typename U2>
struct SubtractInterface<LHS&, RHS,
                    Concepts::DenseMatrix<S1, Orient1, U1>,
                    Concepts::DenseMatrix<S2, Orient2, U2>>
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      subtract(x, swap_sort_order(y));
   }
};

template <typename LHS, typename RHS, typename S1, typename Orient1, typename U1,
          typename S2, typename Orient2, typename U2>
struct SubtractInterface<LHS&, RHS,
                    Concepts::DenseMatrix<S1, Orient1, U1>,
                    Concepts::CompressedOuterMatrix<S2, Orient2, U2>>
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      subtract(swap_sort_order(x), y);
   }
};

// same-major

template <typename LHS, typename RHS, typename Orient,
          typename S1, typename U1, typename S2, typename U2>
struct SubtractInterface<LHS&, RHS,
                    Concepts::DenseMatrix<S1, Orient, U1>,
                    Concepts::DenseMatrix<S2, Orient, U2>>
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      PRECONDITION_EQUAL(size1(x), size1(y));
      PRECONDITION_EQUAL(size2(x), size2(y));
      iter_subtract(iterate(x), const_iterate(y));
   }
};

template <typename LHS, typename RHS, typename Orient,
          typename S1, typename U1, typename S2, typename U2>
struct SubtractInterface<LHS&, RHS,
                    Concepts::CompressedOuterMatrix<S1, Orient, U1>,
                    Concepts::CompressedOuterMatrix<S2, Orient, U2>>
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      PRECONDITION_EQUAL(size1(x), size1(y));
      PRECONDITION_EQUAL(size2(x), size2(y));
      typename Iterate<RHS>::result_type r = iterate(y);
      while (r)
      {
         subtract_element(vector_view(x), r.index(), *r);
         ++r;
      }
   }
};

template <typename LHS, typename RHS, typename Orient,
          typename S1, typename U1, typename S2, typename U2>
struct SubtractInterface<LHS&, RHS,
                    Concepts::DenseMatrix<S1, Orient, U1>,
                    Concepts::CompressedOuterMatrix<S2, Orient, U2>>
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      PRECONDITION_EQUAL(size1(x), size1(y));
      PRECONDITION_EQUAL(size2(x), size2(y));
      typedef typename const_iterator<RHS>::type outer_iterator;
      typedef typename const_iterator<outer_iterator>::type inner_iterator;
      outer_iterator r = iterate(y);
      while (r)
      {
         inner_iterator s = iterate(r);
         while (s)
         {
            subtract_element(x, s.index1(), s.index2(), *s);
            ++s;
         }
         ++r;
      }
   }
};

// multiply

// multiply

template <typename LHS, typename RHS, typename S1, typename U1, typename RHSi>
struct MultiplyInterface<LHS&, RHS, Concepts::LocalMatrix<S1, U1>, AnyScalar<RHSi> >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef RHS const& second_argument_type;
   void operator()(LHS& x, RHS const& y) const
   {
      typename Iterate<LHS&>::result_type i = iterate(x);
      while (i)
      {
         typename inner_iterator<LHS&>::type j = iterate(i);
         while (j)
         {
            multiply(*j, y);
            ++j;
         }
         ++i;
      }
   }
};//
// stream output
//

template <typename T, typename S, typename U>
struct StreamInsert<T, Concepts::SparseMatrix<S, U>>
{
   typedef std::ostream& result_type;
   typedef std::ostream& first_argument_type;
   typedef T second_argument_type;
   result_type operator()(std::ostream& out, T const& x) const;
};

template <typename T, typename S, typename U>
std::ostream&
StreamInsert<T, Concepts::SparseMatrix<S, U>>::operator()(std::ostream& out, T const& x) const
{
   typename const_iterator<T>::type I = iterate(x);
   if (!I)
   {
      return out << "((empty))";
   }
   // else
   out << '(';
   bool first = true;
   while (I)
   {
      typename const_inner_iterator<T>::type Inner = iterate(I);
      while (Inner)
      {
         if (!first)
         {
            out << ", ";
         }
         else first = false;
         out << "(index=(" << Inner.index1() << "," << Inner.index2()
             << "), value=" << *Inner << ")";
         ++Inner;
      }
      ++I;
   }
   out << ')';
   return out;
}

template <typename T, typename S, typename U>
struct StreamInsert<T, Concepts::DenseRowMajorMatrix<S, U>>
{
   typedef std::ostream& result_type;
   typedef std::ostream& first_argument_type;
   typedef T second_argument_type;
   result_type operator()(std::ostream& out, T const& x) const;
};

template <typename T, typename S, typename U>
std::ostream&
StreamInsert<T, Concepts::DenseRowMajorMatrix<S, U>>::operator()(std::ostream& out, T const& x) const
{
   typename const_iterator<T>::type I = iterate(x);
   if (!I)
   {
      return out << "((empty))";
   }
   // else
   out << '(';
   bool first = true;
   while (I)
   {
      if (!first) out << '\n';
      else first = false;
      out << *I;
      ++I;
   }
   out << ')';
   return out;
}

template <typename T, typename S, typename U>
struct StreamInsert<T, Concepts::CompressedRowMajorMatrix<S, U>>
{
   typedef std::ostream& result_type;
   typedef std::ostream& first_argument_type;
   typedef T second_argument_type;
   result_type operator()(std::ostream& out, T const& x) const;
};

template <typename T, typename S, typename U>
std::ostream&
StreamInsert<T, Concepts::CompressedRowMajorMatrix<S, U>>::operator()(std::ostream& out, T const& x) const
{
   typename const_iterator<T>::type I = iterate(x);
   if (!I)
   {
      return out << "((empty))";
   }
   // else
   out << '(';
   bool first = true;
   while (I)
   {
      if (!first) out << '\n';
      else first = false;
      out << "row=" << I.index() << ", value=" << *I;
      ++I;
   }
   out << ')';
   return out;
}


template <typename T, typename S, typename U>
struct StreamInsert<T, Concepts::CompressedColMajorMatrix<S, U>>
{
   typedef std::ostream& result_type;
   typedef std::ostream& first_argument_type;
   typedef T second_argument_type;
   result_type operator()(std::ostream& out, T const& x) const;
};

template <typename T, typename S, typename U>
std::ostream&
StreamInsert<T, Concepts::CompressedColMajorMatrix<S, U>>::operator()(std::ostream& out, T const& x) const
{
   typename const_iterator<T>::type I = iterate(x);
   if (!I)
   {
      return out << "((empty))";
   }
   // else
   out << '(';
   bool first = true;
   while (I)
   {
      if (!first) out << '\n';
      else first = false;
      out << "col=" << I.index() << ", value=" << *I;
      ++I;
   }
   out << ')';
   return out;
}

template <typename T, typename S, typename U>
struct StreamInsert<T, Concepts::DenseColMajorMatrix<S, U>>
{
   typedef std::ostream& result_type;
   typedef std::ostream& first_argument_type;
   typedef T second_argument_type;
   result_type operator()(std::ostream& out, T const& x) const
   {
      return out << swap_sort_order(x);
   }
};

template <typename T, typename S, typename U>
struct StreamInsert<T, Concepts::MatrixExpression<S,U>>
{
   typedef std::ostream& result_type;
   typedef std::ostream& first_argument_type;
   typedef T second_argument_type;
   result_type operator()(std::ostream& out, T const& x) const
   {
      return out << eval_expression(x);
   }
};

} // namespace LinearAlgebra

#endif
