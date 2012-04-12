// -*- C++ -*- $Id$

#if !defined(MATRIX_UTILITY_H_JSHDCUIH78943QY7YPO7YP89)
#define MATRIX_UTILITY_H_JSHDCUIH78943QY7YPO7YP89

#include "matrix.h"
#include <stdlib.h>
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
  return double(rand()) / RAND_MAX;
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
   double r1 = random<double>();
   double r2 = random<double>();
   double u = -std::log(1.0 - r1);
   double rho = std::sqrt(2.0 * u);
   double theta =2*math_const::pi * r2;
   double x = rho*std::cos(theta);           // y = rho*sin(theta), not a bad idea to discard anyway
   return x;
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
