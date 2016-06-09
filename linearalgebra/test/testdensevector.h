// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/test/testdensevector.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "testvectorgeneric.h"

template <typename VecType>
void test_dense_vector()
{
   test_ctor<VecType>();

   typename LinearAlgebra::interface<VecType>::value_type x = 2;

   test_dense_single<VecType>(x);

   VecType y(10,x);
   test_assign(y);

   test_double_negate<VecType>();
   test_real_scalar_vector_nop<VecType>(y);
   test_real_dense_vector<VecType>();
}

template <typename VecType>
void test_real_vector()
{
   VecType v1(10);
   VecType v2(10);

   assign(v1, LinearAlgebra::Range(0,10));
   assign(v2, LinearAlgebra::Range(20,30));
   for (int i = 0; i < 10; ++i)
   {
      CHECK_EQUAL(get_element(v1,i), i);
      CHECK_EQUAL(get_element(v2,i), i+20);
   }

   test_addition(v1,v2);
   test_subtraction(v1,v2);
   test_real_abs(v1);
   test_real_abs(-v1);
   test_equal_to(v1);
   test_minmax(-v1);
}
