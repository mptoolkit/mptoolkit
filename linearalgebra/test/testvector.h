// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/test/testvector.h
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

template <typename VecType>
void test_vector()
{
   VecType v1;
   CHECK_EQUAL(v1.size(), 0);
   CHECK_EQUAL(size(v1), 0);
   CHECK_EQUAL(nnz(v1), 0);

   VecType v2(3);
   CHECK_EQUAL(v2.size(), 3);
   CHECK_EQUAL(size(v2), 3);
   CHECK_EQUAL(nnz(v2), 0);

   VecType v3(4);
   set_element(v3, 0, 0.7);
   CHECK_EQUAL(nnz(v3), 1);
   add_element(v3, 1, 1);
   sub_element(v3, 2, -0.7);
   CHECK_EQUAL(nnz(v3), 3);

   VecType const& v3c(v3);  // to call const members
   CHECK_EQUAL(v3c[0], 0.7);
   CHECK_EQUAL(v3c[1], 1);
   CHECK_EQUAL(v3c[2], 0.7);
   CHECK_EQUAL(v3c[3], 0);    // will not insert a zero element
   CHECK_EQUAL(nnz(v3), 3);

   CHECK_EQUAL(v3[0], 0.7);
   CHECK_EQUAL(v3[1], 1);
   CHECK_EQUAL(v3[2], 0.7);
   CHECK_EQUAL(v3[3], 0);       // non-const; will insert a zero element
   CHECK_EQUAL(nnz(v3), 4);
   v3[3] = 5;
   CHECK_EQUAL(v3[3], 5);

   zero_element(v3, 0);
   CHECK_EQUAL(size(v3), 4);
   CHECK_EQUAL(nnz(v3), 3);
   CHECK_EQUAL(v3c[0], 0);

   VecType v4(-v3);
   CHECK_EQUAL(size(v4), 4);
   CHECK_EQUAL(nnz(v4), 3);
   CHECK_EQUAL(v4[1], -1);
   CHECK_EQUAL(v4[2], -0.7);
   CHECK_EQUAL(v4[3], -5);

   // check that negating twice gets us back to the same type,
   // so the double negation is absorbed.
   CHECK(typeid(transform(transform(v3, LinearAlgebra::Negate<double>()),
                          LinearAlgebra::Negate<double>()))
         == typeid(v3));

   // operator== for rhs proxy
   CHECK_EQUAL(v4, -v3);

   // sanity check
   CHECK(v4 != v3)(v4)(v3);

   // real, conj, transpose, herm should all be identity transformations
   // for real vectors
   CHECK_EQUAL(real(v3), -v4);

   CHECK(typeid(real(v3)) == typeid(v3));
   CHECK(typeid(conj(v3)) == typeid(v3));
   CHECK(typeid(transpose(v3)) == typeid(v3));
   CHECK(typeid(herm(v3)) == typeid(v3));

   // in addition, they should be the same object
   CHECK(&real(v3) == &v3)(&real(v3))(&v3);
   CHECK(&conj(v3) == &v3)(&conj(v3))(&v3);
   CHECK(&transpose(v3) == &v3)(&transpose(v3))(&v3);
   CHECK(&herm(v3) == &v3)(&herm(v3))(&v3);

   // imag part should be zero
   CHECK_EQUAL(norm_1(imag(v3)), 0);

   // check assignment

   v3 = v4;
   CHECK_EQUAL(size(v3), 4);
   CHECK_EQUAL(nnz(v3), 3);
   CHECK_EQUAL(v3[1], -1);
   CHECK_EQUAL(v3[2], -0.7);
   CHECK_EQUAL(v3[3], -5);

   assign(v3, -real(v4) * 2);
   CHECK_EQUAL(size(v3), 4);
   CHECK_EQUAL(nnz(v3), 3);
   CHECK_EQUAL(v3[1], 2);
   CHECK_EQUAL(v3[2], 1.4);
   CHECK_EQUAL(v3[3], 10);

   add_element(v3, 0, 1);
   set_element(v3, 1, -1.4);
   sub_element(v3, 2, 2.4);
   CHECK_EQUAL(nnz(v3), 4);
   CHECK_EQUAL(v3[0], 1);
   CHECK_EQUAL(v3[1], -1.4);
   CHECK_EQUAL(v3[2], -1);
   CHECK_EQUAL(v3[3], 10);

   TRACE(inner_prod(v3, v3));
   TRACE(v3)(v4)(VecType(v3+v4));
   TRACE(-v4);
   TRACE(v4)(VecType(v3-v4));
}
