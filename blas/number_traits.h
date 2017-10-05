// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// blas/matrix.h
//
// Copyright (C) 2017 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

// number_traits schematic
//
// The BLAS trans parameter follows the usual BLAS conventions
// 'N' - normal
// 'T' - transpose
// 'C' - Hermitian conjugate
// 'R' - complex conjugate (BLAS extension)
//
// number_traits<T>::blas_XXXX(char c) returns the trans parameter for
// the given transformation, XXXX = trans, conj, herm.
//
// template <typename NumberType>
// struct number_traits
// {
//    using type = NumberType;
//    static constexpr type zero();
//    static constexpr type identity();
//    static blas_trans(char c);
//    static blas_conj(char c);
//    static blas_herm(char c);
// };

#if !defined(MPTOOLKIT_BLAS_NUMBER_TRAITS_H)
#define MPTOOLKIT_BLAS_NUMBER_TRAITS_H

namespace blas
{

namespace detail
{

struct blas_trans_real
{
   static char blas_trans(char c)
   {
      DEBUG_CHECK(c == 'N' || c == 'T');
      return (c == 'N' ? 'T' : 'N');
   }

   static char blas_conj(char c)
   {
      DEBUG_CHECK(c == 'N' || c == 'T');
      return c;
   }

   static char blas_herm(char c)
   {
      DEBUG_CHECK(c == 'N' || c == 'T');
      return (c == 'N' ? 'T' : 'N');
   }
};

struct blas_trans_complex
{
   static char blas_trans(char c)
   {
      DEBUG_CHECK(c == 'N' || c == 'T' || c == 'C' || c == 'R');
      switch (c)
      {
      case 'N' : return 'T';
      case 'T' : return 'N';
      case 'C' : return 'R';
      case 'R' : return 'C';
      }
   }

   static char blas_conj(char c)
   {
      DEBUG_CHECK(c == 'N' || c == 'T');
      switch (c)
      {
      case 'N' : return 'R';
      case 'T' : return 'C';
      case 'C' : return 'T';
      case 'R' : return 'N';
      }
   }

   static char blas_herm(char c)
   {
      DEBUG_CHECK(c == 'N' || c == 'T');
      switch (c)
      {
      case 'N' : return 'C';
      case 'T' : return 'R';
      case 'C' : return 'N';
      case 'R' : return 'T';
      }
   }
};

} // namespace detail

//
// number_traits
// basic information about scalar types that can be put in dense matrices
//

template <typename T>
struct number_traits : public detail::blas_trans_real
{
   using type = T;

   static constexpr T zero() { return 0; }
   static constexpr T identity() { return 1; }
};

template <>
struct number_traits<float> : public detail::blas_trans_real
{
   using type = float;

   static constexpr float zero() { return 0.0f; }
   static constexpr float identity() { return 1.0f; }
};

template <>
struct number_traits<double> : public detail::blas_trans_real
{
   using type = double;

   static constexpr double zero() { return 0.0; }
   static constexpr double identity() { return 1.0; }
};

template <>
struct number_traits<std::complex<float>> : public detail::blas_trans_complex
{
   using type = std::complex<float>;

   static constexpr std::complex<float> zero() { return {0.0f,0.0f}; }
   static constexpr std::complex<float> identity() { return {1.0f,0.0f}; }
};

template <>
struct number_traits<std::complex<double>> : public detail::blas_trans_complex
{
   using type = std::complex<double>;

   static constexpr std::complex<double> zero() { return {0.0,0.0}; }
   static constexpr std::complex<double> identity() { return {1.0,0.0}; }
};

#if defined(HAVE_FLOAT128)
template <>
struct number_traits<float128> : public blas_trans_real
{
   using type = float128;

   static constexpr float128 zero() { return 0.0Q; }
   static constexpr float128 identity() { return 1.0Q; }
};

template <>
struct number_traits<std::complex<float128>> : public detail::blas_trans_complex
{
   using type = std::complex<float128>;

   static constexpr std::complex<float128> zero() { return {0.0Q,0.0Q}; }
   static constexpr std::complex<float128> identity() { return {1.0Q,0.0Q}; }
};
#endif

} // namespace blas

#endif
