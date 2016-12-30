// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/associate.h
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
  some experimental stuff on how to implement associativity of operators
*/

//
// AssociateRight
//
// a traits type that determines whether (x op1 (y op2 z)) can be transformed into
// (x op1' y) op2' z.
// Whether this *should* be done in a given expression, is not determined here,
// except transformations that could be done but would result in assymptotically worse
// expressions are (or should) not be done.  For example, matrix * (matrix * vector) should
// never be transformed into (matrix * matrix) * vector.
//

template <typename BinaryFunctor1, typename BinaryFunctor2, typename Enable = void>
struct AssociateRight : public boost::mpl::false_
{
};

template <typename BinaryFunctor1, typename BinaryFunctor2, typename Enable = void>
struct AssociateLeft : public boost::mpl::false_
{
};

// multiplication is associative if there is a promote traits
template <typename T1, typename T23, typename T2, typename T3>
struct AssociateRight<BinaryOperator<Multiplication, T1, T23>, BinaryOperator<Multiplication, T2, T3>,
                      boost::enable_if<boost::and_<has_value_type<PromoteTraits<T1, T23> >,
                                                   has_value_type<PromoteTraits<T2, T3> > > >::type>
  : public boost::mpl
{
   typedef BinaryOpertor<Multiplication, T1, T2> LeftOperator;
   typedef BinaryOperator<Multiplication, typename LeftOperator::value_type, T3> RightOperator;
};

// multiplication is associative if there is a promote traits
template <typename T1, typename T2, typename T12, typename T3>
struct AssociateLeft<BinaryOperator<Multiplication, T1, T2>, BinaryOperator<Multiplication, T12, T3>,
                      boost::enable_if<boost::and_<has_value_type<PromoteTraits<T1, T2> >,
                                                   has_value_type<PromoteTraits<T12, T3> > > >::type>
  : public boost::mpl
{
   typedef BinaryOperator<Multiplication, T2, T3> RightOperator;
   typedef BinaryOperator<Multiplication, T1, typename RightOperator::value_type> LeftOperator;
};





// transform scalar * (scalar * matrix) into (scalar * scalar) * matrix
template <typename T1, typename T23, typename T2, typename T3>
struct AssociateRight<ScalarMatrixProduct<T1, T23>, ScalarMatrixProduct<T2, T3>,
 boost::enable_if<AssociateRight<BinaryOperator<Multiplication, T1, T2>,
                                 BinaryOperator<Multiplication,
                                                typename BinaryOperator<Multiplication, T1, T2>::value_type,
                                                T3> > >::type>
  : public boost::mpl::true_
{
   typedef BinaryOperator<Multiplication, T1, T2> LeftOperator;
   typedef ScalarMatrixProduct<typename LeftOperator::value_type, T3> RightOperator;
};












template <typename Tag, typename T1, typename T2, typename T3>
struct AssociativeTraits<T1, Multiplication, T2, Multiplication, T3,
                         typename boost::enable_if<
                            boost::and_<
                               has_value_type<PromoteTraits<T1, T2> >,
                               has_value_type<PromoteTraits<T2, T3> >,
                            >
                         >::type>
{
   static bool const value = true;
};



// multiply by scalar associate iff the underlying types do
template <typename T1, typename T2, typename T3>
struct AssociativityTraits<T1, Multiplication, T2, ScalarMatrixMultiplication, T3>
  : public AssociativityTraits<T1, Multiplication, T2, Multiplication, typename T3::value_type>
{
};

template <typename T1, typename T2, typename T3>
struct AssociativityTraits<T1, MatrixScalarMultiplication, T2, Multiplication, T3>
  : public AssociativityTraits<typename T1::value_type, Multiplication, T2, Multiplication, T3>
{
};

// matrix-matrix multiply associates iff the underlying types do
template <typename T1, typename T2, typename T3>
struct AssociativityTraits<T1, MatrixMatrixMultiplication, T2, MatrixMatrixMultiplication, T3>
  : public AssociativityTraits<typename T1::value_type, Multiplication,
                               typename T2::value_type, Multiplication,
                               typename T3::value_type>
{
};

// matrix-matrix multiply associates iff the underlying types do
template <typename T1, typename T2, typename T3>
struct AssociativityTraits<T1, MatrixMatrixMultiplication, T2, MatrixMatrixMultiplication, T3>
  : public AssociativityTraits<typename T1::value_type, Multiplication,
                               typename T2::value_type, Multiplication,
                               typename T3::value_type>
{
};

template <typename T1, typename T2, typename T3>
struct AssociativityTraits<T1, MatrixMatrixMultiplication, T2, MatrixScalarMultiplication, T3>
  : public AssociativityTraits<typename T1::value_type, Multiplication,
                               typename T2::value_type, Multiplication,
                               T3>
{
};

template <typename T1, typename T2, typename T3>
struct AssociativityTraits<T1, ScalarMatrixMultiplication, T2, MatrixMatrixMultiplication, T3>
  : public AssociativityTraits<T1, Multiplication,
                               typename T2::value_type, Multiplication,
                               typename T3::value_type>
{
};
