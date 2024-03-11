// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/scaledmatrix.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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

#if !defined(SCALEDMATRIX_H_SDCHJKSDHUIRYTYT78YHIULRELE)
#define SCALEDMATRIX_H_SDCHJKSDHUIRYTYT78YHIULRELE

#include "identitymatrix.h"
#include "directproduct.h"

namespace LinearAlgebra
{

template <typename Mat, typename Scalar = typename Mat::value_type>
class ScaledMatrix : public GenericMatrix<typename BinaryOperator<Multiplication,
                                                                  Scalar,
                                                                  typename Mat::value_type>::value_type,
                                          ScaledMatrix<Mat, Scalar> >
{
   public:
      typedef Mat base_type;
      typedef typename BinaryOperator<Multiplication,
                                      Scalar,
                                      typename Mat::value_type>::value_type value_type;

      typedef Scalar scalar_type;

      ScaledMatrix(base_type const& b, scalar_type const& s);

      base_type& base() { return base_; }
      base_type const& base() const { return base_; }

      scalar_type& scalar() { return scalar_; }
      scalar_type const& scalar() const { return scalar_; }

   private:
       base_type base_;
       scalar_type scalar_;
};

template <typename Mat, typename Scalar>
struct MatrixTraits<ScaledMatrix<Mat, Scalar> >
{
   typedef ScaledMatrix<Mat, Scalar> expression_type;
};


template <typename B, typename S, typename T>
struct MatrixMatrixProduct<ScaledMatrix<IdentityMatrix<B>, S>,
                        ScaledMatrix<IdentityMatrix<B>, T> >
{
   typedef typename BinaryOperator<Multiplication, S, T>::value_type ValType;
   typedef ScaledMatrix<IdentityMatrix<B>, ValType> result_type;
   typedef result_type value_type;

   typedef ScaledMatrix<IdentityMatrix<B>, S> first_argument_type;
   typedef ScaledMatrix<IdentityMatrix<B>, T> second_argument_type;

   result_type operator()(first_argument_type const& x, second_argument_type const& y)
   {
      DEBUG_PRECONDITION(x.size2() == y.size1());
      return result_type(IdentityMatrix<B>(x.size1(), y.size2()), x.scalar() * y.scalar());
   }
};

template <typename B, typename S, typename T>
struct MatrixScalarProduct<ScaledMatrix<IdentityMatrix<B>, S>, T>
{
   typedef typename BinaryOperator<Multiplication, S, T>::value_type ValType;
   typedef ScaledMatrix<IdentityMatrix<B>, ValType> result_type;
   typedef result_type value_type;

   typedef ScaledMatrix<IdentityMatrix<B>, S> first_argument_type;
   typedef T second_argument_type;

   result_type operator()(first_argument_type const& x, second_argument_type const& y)
   {
      return result_type(x.base(), x.scalar() * y);
   }
};

template <typename B, typename S, typename T>
struct ScalarMatrixProduct<S, ScaledMatrix<IdentityMatrix<B>, T> >
{
   typedef typename BinaryOperator<Multiplication, S, T>::value_type ValType;
   typedef ScaledMatrix<IdentityMatrix<B>, ValType> result_type;
   typedef result_type value_type;

   typedef S first_argument_type;
   typedef ScaledMatrix<IdentityMatrix<B>, T> second_argument_type;

   result_type operator()(first_argument_type const& x, second_argument_type const& y)
   {
      return result_type(y.base(), x * y.scalar());
   }
};

//
// DirectProduct
//

template <typename B, typename S, typename T>
struct DirectProduct<ScaledMatrix<IdentityMatrix<B>, S>, ScaledMatrix<IdentityMatrix<B>, T> >
{
   typedef typename BinaryOperator<Multiplication, S, T>::value_type ValType;
   typedef ScaledMatrix<IdentityMatrix<B>, ValType> result_type;
   typedef result_type value_type;

   typedef ScaledMatrix<IdentityMatrix<B>, S> first_argument_type;
   typedef ScaledMatrix<IdentityMatrix<B>, T> second_argument_type;

   result_type operator()(first_argument_type const& x, second_argument_type const& y)
   {
      return result_type(IdentityMatrix<B>(x.size1() * y.size1()), x.scalar() * y.scalar());
   }
};

template <typename B, typename S, typename T>
struct DirectProduct<S, ScaledMatrix<IdentityMatrix<B>, T>,
                     typename boost::enable_if<is_matrix<S> >::type>
{
   typedef S first_argument_type;
   typedef typename S::value_type SValType;
   typedef ScaledMatrix<IdentityMatrix<B>, T> second_argument_type;

   typedef typename BinaryOperator<Multiplication, SValType, T>::value_type ValType;
   typedef Matrix<ValType> value_type;
   typedef value_type result_type;

   result_type operator()(first_argument_type const& x, second_argument_type const& y) const
   {
      size_type m2s1 = y.size1();
      size_type m2s2 = y.size2();
      value_type Res(x.size1() * y.size1(), x.size2() * y.size2(), 0);
      for (int i = 0; i < x.size1(); ++i)
      {
         for (int j = 0; j < x.size2(); ++j)
         {
            Res.range1(i * m2s1, (i+1) * m2s1).range2(j * m2s2, (j+1) * m2s2)
               = x(i,j) * y;
         }
      }
      return Res;
   }
};

// FIXME: hmm, the scalar should appear on the other side of the multiply here
template <typename B, typename S>
struct Conj<ScaledMatrix<B, S> >
{
   typedef ScaledMatrix<B, S> argument_type;
   typedef typename Conj<B>::value_type ValType;
   typedef typename Conj<S>::value_type SValType;

   typedef ScaledMatrix<ValType, SValType> result_type;
   typedef result_type value_type;

   result_type operator()(argument_type const& x) const
   {
      return result_type(conj(x.base()), conj(x.scalar()));
   }
};

template <typename B, typename S, typename T>
struct DirectProduct<ScaledMatrix<IdentityMatrix<B>, S>, T,
                     typename boost::enable_if<is_matrix<T> >::type>
{
   typedef ScaledMatrix<IdentityMatrix<B>, S> first_argument_type;
   typedef T second_argument_type;
   typedef typename T::value_type TValType;

   typedef typename BinaryOperator<Multiplication, S, TValType>::value_type ValType;
   typedef Matrix<ValType> value_type;
   typedef value_type result_type;

   result_type operator()(first_argument_type const& x, second_argument_type const& y) const
   {
      size_type m2s1 = y.size1();
      size_type m2s2 = y.size2();
      value_type Res(x.size1() * m2s1, x.size2() * m2s2, 0);
      for (int i = 0; i < x.size1(); ++i)
      {
         Res.range1(i * m2s1, (i+1) * m2s1).range2(i * m2s2, (i+1) * m2s2)
            = x.scalar() * y;
      }
      return Res;
   }
};

template <class Scalar, class Derived>
template <typename B, typename T>
void MatrixExpression<Scalar, Derived>::assign(ScaledMatrix<IdentityMatrix<B>, T> const& V)
{
   CHECK_EQUAL(this->size1(), V.size1());
   CHECK_EQUAL(this->size2(), V.size2());

   this->fill(T(B()));
   int i = 0;
   for (iterator1 I = this->begin(); I != this->end(); ++I, ++i)
   {
      I[i] = V.scalar();
   }
}

template <class Scalar, class Derived>
template <typename B, typename T>
void MatrixExpression<Scalar, Derived>::add(ScaledMatrix<IdentityMatrix<B>, T> const& V)
{
   CHECK_EQUAL(this->size1(), V.size1());
   CHECK_EQUAL(this->size2(), V.size2());

   int i = 0;
   for (iterator1 I = this->begin(); I != this->end(); ++I, ++i)
   {
      I[i] += V.scalar();
   }
}

template <class Scalar, class Derived>
template <typename B, typename T>
void MatrixExpression<Scalar, Derived>::subtract(ScaledMatrix<IdentityMatrix<B>, T> const& V)
{
   CHECK_EQUAL(this->size1(), V.size1());
   CHECK_EQUAL(this->size2(), V.size2());

   int i = 0;
   for (iterator1 I = this->begin(); I != this->end(); ++I, ++i)
   {
      I[i] += V.scalar();
   }
}


} // namespace LinearAlgebra

#endif
