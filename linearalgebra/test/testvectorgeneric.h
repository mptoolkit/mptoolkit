// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/test/testvectorgeneric.h
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

template <typename T>
void test_ctor()
{
   // default ctor
   T x;
   CHECK_EQUAL(size(x), 0);

   // copy ctor
   T y(x);
   CHECK_EQUAL(size(y), 0);

   // assignment
   x = y;
   CHECK_EQUAL(size(x), 0);
}

template <typename T, typename Value>
void test_dense_single(Value const& x)
{
   // test size 1 vectors
   T v(1);
   CHECK_EQUAL(size(v), 1);

   set_element(v,0,x);
   CHECK_EQUAL(get_element(v,0), x);

   T w(1,x);
   CHECK_EQUAL(size(w), 1);
   CHECK_EQUAL(get_element(v,0), get_element(w,0));
   CHECK_EQUAL(v,w);

   T z(1);
   fill(z, x);
   CHECK_EQUAL(size(z), 1);
   CHECK_EQUAL(get_element(z,0), x);
}

template <typename T>
void test_assign(T const& x)
{
   T y;
   CHECK_EQUAL(size(y), 0);

   y = x;
   CHECK_EQUAL(size(y), size(x));
   CHECK_EQUAL(y,x);

   TRACE(size(x));
   T z(size(x));
   assign(z, x);
   CHECK_EQUAL(size(z), size(x));
   CHECK_EQUAL(z,x);
}

// check that double negation converts to a no-op
template <typename T>
void test_double_negate()
{
   T x;
   CHECK(typeid(transform(transform(x, LinearAlgebra::Negate<double>()),
                                LinearAlgebra::Negate<double>())) ==
         typeid(x));
}

template <typename T>
void test_real_scalar_vector_nop(T const& v)
{
   CHECK(typeid(real(v)) == typeid(v));
   CHECK(typeid(conj(v)) == typeid(v));
   CHECK(typeid(transpose(v)) == typeid(v));
   CHECK(typeid(herm(v)) == typeid(v));

   // imag part should be zero
   CHECK_EQUAL(norm_1(imag(v)), 0);

   // in addition, they should be the same object
   CHECK_EQUAL(&real(v), &v);
   CHECK_EQUAL(&conj(v), &v);
   CHECK_EQUAL(&transpose(v), &v);
   CHECK_EQUAL(&herm(v), &v);
}

template <typename VecType>
void test_real_dense_vector()
{
   VecType v1;
   CHECK_EQUAL(size(v1), 0);

   VecType v2(3);
   CHECK_EQUAL(size(v2), 3);

   VecType v3(3, 0.7);
   CHECK_EQUAL(size(v3), 3);
   CHECK_EQUAL(v3[0], 0.7);
   CHECK_EQUAL(v3[1], 0.7);
   CHECK_EQUAL(v3[2], 0.7);

   VecType v4(size(v3));
   assign(v4, -v3);
   CHECK_EQUAL(size(v4), 3);
   CHECK_EQUAL(v4[0], -0.7);
   CHECK_EQUAL(v4[1], -0.7);
   CHECK_EQUAL(v4[2], -0.7);

   // operator== for rhs proxy
   CHECK_EQUAL(v4, -v3);

   // sanity check
   CHECK(v4 != v3)(v4)(v3);

   // check assignment
   v3 = v4;
   CHECK_EQUAL(size(v3), 3);
   CHECK_EQUAL(v3[0], -0.7);
   CHECK_EQUAL(v3[1], -0.7);
   CHECK_EQUAL(v3[2], -0.7);

   assign(v3, -real(v4) * 2);
   CHECK_EQUAL(size(v3), 3);
   CHECK_EQUAL(v3[0], 1.4);
   CHECK_EQUAL(v3[1], 1.4);
   CHECK_EQUAL(v3[2], 1.4);

   // nnz
   CHECK_EQUAL(nnz(v3), 3);
   add_element(v3, 0, 1);
   set_element(v3, 1, -1.4);
   subtract_element(v3, 2, 2.4);
   CHECK_EQUAL(v3[0], 2.4);
   CHECK_EQUAL(v3[1], -1.4);
   CHECK_EQUAL(v3[2], -1);
}

template <typename T>
void test_addition(T const& x, T const& y)
{
   T z(size(x));
   assign(z,x+y);

   for (int i = 0; i < int(size(x)); ++i)
   {
      CHECK_EQUAL(get_element(z,i), get_element(x,i) + get_element(y,i));
   }
}

template <typename T>
void test_subtraction(T const& x, T const& y)
{
   T z(size(x));
   assign(z,x-y);

   for (int i = 0; i < int(size(x)); ++i)
   {
      CHECK_EQUAL(get_element(z,i), get_element(x,i) - get_element(y,i));
   }
}

template <typename T>
void test_minmax(T const& x)
{
   typename LinearAlgebra::make_value<typename
      LinearAlgebra::GetVectorElement<T>::result_type>::type m = min(x);
   typename LinearAlgebra::make_value<typename
      LinearAlgebra::GetVectorElement<T>::result_type>::type M = max(x);

   CHECK_EQUAL(min(x), -max(-x));
   CHECK_EQUAL(max(x), -min(-x));

   typename LinearAlgebra::const_iterator<T>::type I = iterate(x);
   while (I)
   {
      CHECK(!(*I < m));
      CHECK(!(M < *I));
      ++I;
   }
}

template <typename T>
void test_real_abs(T const& x)
{
   typedef typename LinearAlgebra::make_value<typename
      LinearAlgebra::Abs<T>::result_type>::type AbsType;

   using LinearAlgebra::equal;
   using LinearAlgebra::abs;

   AbsType v = abs(x);
   CHECK_CLOSE(norm_1(x), norm_1(v));

   typename LinearAlgebra::interface<AbsType>::value_type r = 0;
   typename LinearAlgebra::const_iterator<AbsType>::type I = iterate(v);
   while (I)
   {
      r += *I;
      ++I;
   }
   CHECK_CLOSE(norm_1(x), r);
}

template <typename T>
void test_equal_to(T const& x)
{
   CHECK_EQUAL(x, x);
   CHECK(LinearAlgebra::equal_to(x,x));
}
