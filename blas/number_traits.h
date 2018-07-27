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
//    using real_type = (type, if it is real, otherwise T, if type is complex<T>)
//    using complex_type = (type, if is complex, complex<type>)
//    static constexpr type zero();
//    static constexpr type identity();
//    static blas_trans(char c);
//    static blas_conj(char c);
//    static blas_herm(char c);
// };

#if !defined(MPTOOLKIT_BLAS_NUMBER_TRAITS_H)
#define MPTOOLKIT_BLAS_NUMBER_TRAITS_H

#include "common/trace.h"
#include <cmath>
#include <type_traits>
#include "common/types.h"

namespace blas
{

struct cpu_tag;

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
      return 'x';
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
      return 'x';
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
      return 'x';
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
   using real_type = T;
   using complex_type = std::complex<T>;

   static constexpr T zero() { return 0; }
   static constexpr T identity() { return 1; }
};

template <>
struct number_traits<float> : public detail::blas_trans_real
{
   using type = float;
   using real_type = float;
   using complex_type = std::complex<float>;

   static constexpr float zero() { return 0.0f; }
   static constexpr float identity() { return 1.0f; }
};

template <>
struct number_traits<double> : public detail::blas_trans_real
{
   using type = double;
   using real_type = double;
   using complex_type = std::complex<double>;

   static constexpr double zero() { return 0.0; }
   static constexpr double identity() { return 1.0; }
};

template <>
struct number_traits<std::complex<float>> : public detail::blas_trans_complex
{
   using type = std::complex<float>;
   using real_type = float;
   using complex_type = std::complex<float>;

   static constexpr std::complex<float> zero() { return {0.0f,0.0f}; }
   static constexpr std::complex<float> identity() { return {1.0f,0.0f}; }
};

template <>
struct number_traits<std::complex<double>> : public detail::blas_trans_complex
{
   using type = std::complex<double>;
   using real_type = double;
   using complex_type = std::complex<double>;

   static constexpr std::complex<double> zero() { return {0.0,0.0}; }
   static constexpr std::complex<double> identity() { return {1.0,0.0}; }
};

#if defined(HAVE_FLOAT128)
template <>
struct number_traits<float128> : public blas_trans_real
{
   using type = float128;
   using real_type = float128;
   using complex_type = std::complex<float128>;

   static constexpr float128 zero() { return 0.0Q; }
   static constexpr float128 identity() { return 1.0Q; }
};

template <>
struct number_traits<std::complex<float128>> : public detail::blas_trans_complex
{
   using type = std::complex<float128>;
   using real_type = float128;
   using complex_type = std::complex<float128>;

   static constexpr std::complex<float128> zero() { return {0.0Q,0.0Q}; }
   static constexpr std::complex<float128> identity() { return {1.0Q,0.0Q}; }
};
#endif

// some miscellaneous functions for scalars

inline
void add_inner_prod(double x, double y, double& r)
{
   r += x*y;
}

inline
void add_inner_prod(std::complex<double> x, std::complex<double> y, std::complex<double>& r)
{
   r += std::conj(x)*y;
}

namespace detail
{

template <typename T, typename Enable = void>
struct remove_proxy_helper
{
   using type = T;
};

template <typename T>
struct void_
{
   using type = void;
};

template <typename T>
struct remove_proxy_helper<T, typename void_<typename T::remove_proxy_t>::type>
{
   using type = typename T::remove_proxy_t;
};

template <typename T, typename Enable = void>
struct tag_of_helper
{
   using type = cpu_tag;
};

template <typename T>
struct tag_of_helper<T, typename void_<typename T::tag_type>::type>
{
   using type = typename T::tag_type;
};

} // namespace detail

//
// tag_of
//
// Metafunction to get the storage tag of an arbitrary type.
//

template <typename T>
struct tag_of : detail::tag_of_helper<T>
{
};

template <typename T>
using tag_of_t = typename tag_of<T>::type;

//
// remove_proxy
//
// Metafunction to get the 'base' type of a proxy reference.

template <typename T>
struct remove_proxy : public detail::remove_proxy_helper<T> {};

template <typename T>
using remove_proxy_t = typename remove_proxy<T>::type;

//
// is_numeric
//
// metafunction to determine if a type is a scalar type (ie, an arithmetic type
// or a user type that emulates an arithmetic type).
//

template <typename T>
struct is_numeric : std::is_arithmetic<T> {};

template <typename T>
struct is_numeric<T const> : is_numeric<T> {};

template <typename T>
struct is_numeric<T volatile> : is_numeric<T> {};

template <typename T>
struct is_numeric<T&> : is_numeric<T> {};

template <typename T>
struct is_numeric<T&&> : is_numeric<T> {};

template <typename T>
constexpr bool is_numeric_v = is_numeric<T>::value;

template <typename T>
struct is_numeric<std::complex<T>> : is_numeric<T> {};

#if defined(HAVE_FLOAT128)
template <>
struct is_numeric<float128> : std::true_type {};
#endif

template <typename T>
constexpr bool is_numeric_v = is_numeric<T>::value;

// is_real, is_complex

template <typename T>
struct is_real : std::is_floating_point<T> {};

template <typename T>
constexpr bool is_real_v = is_real<T>::value;

template <typename T>
struct is_complex : std::false_type {};

template <typename T>
struct is_complex<std::complex<T>> : std::true_type {};

template <typename T>
constexpr bool is_complex_v = is_complex<T>::value;

// metafunction to get the underlying numeric type (ie, real or complex, of some precision)
// of a matrix or vector

template <typename T>
struct numeric_type_of;

template <typename T>
using numeric_type_of_t = typename numeric_type_of<T>::type;

} // namespace blas

//
// remainder of header is in the global namespace
//

// copy() function for deep copy.  The two generic versions
// do no-operations for:
// (1) the case where the type has a copy constructor
// (2) the case where we are copying a temporary (or moved-from value) and we have a move constructor

// global copy functions
template <typename T>
inline
std::enable_if_t<std::is_copy_constructible<T>::value, T>
copy(T const& x)
{
   return x;
}

// for a temporary, the version that takes a universal reference is found first
// in overload resolution and we can avoid the copy im this case
template <typename T>
inline
std::enable_if_t<std::is_move_constructible<T>::value && !std::is_lvalue_reference<T>::value, T>
copy(T&& x) // x is a universal reference here, which isn't really what we want
{
   return std::move(x);
}

// uninitialized_default_construct doesn't exist intil C++17 so define it here as an extension

namespace stdext
{

template<class T>
void destroy_at(T* p)
{
    p->~T();
}

template<class ForwardIt, class Size>
ForwardIt destroy_n(ForwardIt first, Size n)
{
  for (; n > 0; (void) ++first, --n)
    stdext::destroy_at(std::addressof(*first));
  return first;
}

template<class ForwardIt, class Size>
ForwardIt uninitialized_default_construct_n(ForwardIt first, Size n)
{
    typedef typename std::iterator_traits<ForwardIt>::value_type Value;
    ForwardIt current = first;
    try {
        for (; n > 0 ; (void) ++current, --n) {
            ::new (static_cast<void*>(std::addressof(*current))) Value;
        }
        return current;
    }  catch (...) {
        for (; first != current; ++first) {
            first->~Value();
        }
        throw;
    }
}

} // namespace stdext

#include "number_traits.icc"

#endif
