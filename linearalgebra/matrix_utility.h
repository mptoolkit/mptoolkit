// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/matrix_utility.h
//
// Copyright (C) 2004-2023 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2022 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#if !defined(MATRIX_UTILITY_H_JSHDCUIH78943QY7YPO7YP89)
#define MATRIX_UTILITY_H_JSHDCUIH78943QY7YPO7YP89

#include "matrix.h"
#include "common/randutil.h"
#include <complex>

namespace LinearAlgebra
{

template <typename Func>
Matrix<typename Func::value_type>
generate_matrix(size_type Size1, size_type Size2, Func f = Func())
{
   Matrix<typename Func::value_type> Result(Size1, Size2);
   for (size_type i = 0; i < Size1; ++i)
   {
      for (size_type j = 0; j < Size2; ++j)
      {
         Result(i,j) = f();
      }
   }
   return Result;
}

template <typename Scalar>
Matrix<Scalar>
generate_matrix(size_type Size1, size_type Size2, Scalar (&f)())
{
   Matrix<typename make_value<Scalar>::type> Result(Size1, Size2);
   for (size_type i = 0; i < Size1; ++i)
   {
      for (size_type j = 0; j < Size2; ++j)
      {
         Result(i,j) = f();
      }
   }
   return Result;
}

template <typename Scalar>
Scalar random();

template <>
inline
double random<double>()
{
  return randutil::rand();
}

template <>
inline
std::complex<double> random<std::complex<double> >()
{
   return std::complex<double>(random<double>(), random<double>());
}

// Uniformly distributed random number
template <typename Scalar>
Scalar nrandom();

template <>
inline
double nrandom<double>()
{
   return randutil::randn();
}

template <>
inline
std::complex<double> nrandom<std::complex<double> >()
{
   return math_const::r_1_sqrt_2 * std::complex<double>(nrandom<double>(), nrandom<double>());
}

template <typename Scalar>
Matrix<Scalar> random_matrix(size_type Size1, size_type Size2);  // prototype

template <typename Scalar>
Matrix<Scalar> nrandom_matrix(size_type Size1, size_type Size2);  // prototype

// work-around for icc that doesn't like the function template
// FIXME: this was probably fixed a decade ago ...
#if defined(INTEL_COMPILER) || defined(__INTEL_COMPILER)

template <>
Matrix<double> random_matrix(size_type Size1, size_type Size2)
{
   return generate_matrix(Size1, Size2, random<double>);
}

template <>
Matrix<std::complex<double> > random_matrix(size_type Size1, size_type Size2)
{
   return generate_matrix(Size1, Size2, random<std::complex<double> >);
}

#else

template <typename Scalar>
Matrix<Scalar> random_matrix(size_type Size1, size_type Size2)
{
   return generate_matrix(Size1, Size2, random<Scalar>);
}

#endif

template <typename Scalar>
Matrix<Scalar> nrandom_matrix(size_type Size1, size_type Size2)
{
   return generate_matrix(Size1, Size2, nrandom<Scalar>);
}

// Random matrix distributed according to the Haar measure.
// For real types, this generates an orthogonal matrix
template <typename Scalar>
Matrix<Scalar>
random_unitary(size_type Size1, size_type Size2)
{
   int sz = std::max(Size1, Size2);
   Matrix<Scalar> X = nrandom_matrix<Scalar>(sz, sz);
   auto M = QR_FactorizeFull(X);
   Matrix<Scalar> Result(Size1, Size2);
   Result = M(LinearAlgebra::range(0,Size1), LinearAlgebra::range(0,Size2));
   return Result;
}

// FIXME: this is a hack in lieu of a proper diagonal matrix type

template <typename T>
Matrix<typename interface<T>::value_type>
diagonal_matrix(T const& v)
{
   size_type sz = size(v);
   Matrix<typename interface<T>::value_type> Result(sz,sz);
   zero_all(Result);
   for (size_type i = 0U; i < sz; ++i)
   {
      Result(i,i) = get_element(v, i);
   }
   return Result;
}

template <typename T>
Vector<T>
diagonal_vector(Matrix<T> const& v)
{
   CHECK_EQUAL(size1(v), size2(v));
   size_type sz = size1(v);
   Vector<T> Result(sz);
   for (size_type i = 0U; i < sz; ++i)
   {
      Result[i] = v(i,i);
   }
   return Result;
}

template <typename T>
Vector<T>
real_diagonal_vector(Matrix<std::complex<T>> const& v)
{
   CHECK_EQUAL(size1(v), size2(v));
   size_type sz = size1(v);
   Vector<T> Result(sz);
   for (size_type i = 0U; i < sz; ++i)
   {
      Result[i] = v(i,i).real();
   }
   return Result;
}

template <typename T>
bool is_symmetric(T const& x, typename boost::enable_if<is_matrix<T> >::type* = 0)
{
   return equal(x, transpose(x));
}

template <typename T>
bool is_hermitian(T const& x, typename boost::enable_if<is_matrix<T> >::type* = 0)
{
   return equal(x, herm(x));
}

template <typename Mat>
void zero_upper_triangular(Mat& M)
{
   for (unsigned i = 0; i < M.size1(); ++i)
   {
      for (unsigned j = i+1; j < M.size2(); ++j)
      {
         M(i,j) = zero<typename interface<Mat>::value_type>();
      }
   }
}

template <typename Mat>
void zero_lower_triangular(Mat& M)
{
   for (unsigned j = 0; j < M.size2(); ++j)
   {
      for (unsigned i = j+1; i < M.size1(); ++i)
      {
         M(i,j) = zero<typename interface<Mat>::value_type>();
      }
   }
}

} // namespace LinearAlgebra

#endif
