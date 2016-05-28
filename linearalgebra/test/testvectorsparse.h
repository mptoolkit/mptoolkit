// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/test/testvectorsparse.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "testvectorgeneric.h"

template <typename VecType>
void test_real_sparse_vector()
{
   test_ctor<VecType>();

   VecType v1;   
   CHECK_EQUAL(size(v1), 0);
   CHECK_EQUAL(nnz(v1), 0);

   VecType v2(3);
   CHECK_EQUAL(size(v2), 3);
   CHECK_EQUAL(nnz(v2), 0);

   VecType v3(4);
   set_element(v3, 0, 0.7);
   CHECK_EQUAL(nnz(v3), 1);
   add_element(v3, 1, 1);
   subtract_element(v3, 2, -0.7);
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

   // operator== for rhs proxy
   CHECK_EQUAL(v4, -v3);

   // sanity check
   CHECK(v4 != v3)(v4)(v3);
 
   test_double_negate<VecType>();
   test_real_scalar_vector_nop(v4);
   
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
   subtract_element(v3, 2, 2.4);
   CHECK_EQUAL(nnz(v3), 4);
   CHECK_EQUAL(v3[0], 1);
   CHECK_EQUAL(v3[1], -1.4);
   CHECK_EQUAL(v3[2], -1);
   CHECK_EQUAL(v3[3], 10);

   // TODO:
   // norm_1, norm_2, norm_2_sq, inner_prod, etc etc
}
