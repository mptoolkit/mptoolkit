// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/test/testdirectsum.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "linearalgebra/vector.h"
#include "linearalgebra/hashvector.h"
#include "linearalgebra/fixedvector.h"

using namespace LinearAlgebra;

int main()
{
   Vector<double> U = range(0,10,2);

   HashVector<double> V(12);
   V[3] = 2;
   V[7] = 5;

   FixedVector<double> W(5, -2);

   Vector<double> Test1 = direct_sum(U,V);
   TRACE(Test1);
   CHECK_EQUAL(size(Test1), 17);
   CHECK_EQUAL(Test1[0], 0);
   CHECK_EQUAL(Test1[1], 2);
   CHECK_EQUAL(Test1[2], 4);
   CHECK_EQUAL(Test1[3], 6);
   CHECK_EQUAL(Test1[4], 8);
   CHECK_EQUAL(Test1[5], 0);
   CHECK_EQUAL(Test1[6], 0);
   CHECK_EQUAL(Test1[7], 0);
   CHECK_EQUAL(Test1[8], 2);
   CHECK_EQUAL(Test1[9], 0);
   CHECK_EQUAL(Test1[10], 0);
   CHECK_EQUAL(Test1[11], 0);
   CHECK_EQUAL(Test1[12], 5);
   CHECK_EQUAL(Test1[13], 0);
   CHECK_EQUAL(Test1[14], 0);
   CHECK_EQUAL(Test1[15], 0);
   CHECK_EQUAL(Test1[16], 0);

   HashVector<double> Test2 = direct_sum(V,V);
   CHECK_EQUAL(size(Test2), 24);
   CHECK_EQUAL(Test2[3], 2);
   CHECK_EQUAL(Test2[7], 5);
   CHECK_EQUAL(Test2[15], 2);
   CHECK_EQUAL(Test2[19], 5);
   CHECK_EQUAL(norm_1(Test2), 14);

   Vector<double> Test3 = direct_sum(U,W);
   CHECK_EQUAL(size(Test3), 10);
   CHECK_EQUAL(norm_1(Test3), 30);
   CHECK_EQUAL(Test3[1], 2);

   Vector<double> Test4 = direct_sum(4 * U, 3 * W, 2 * V);
   CHECK_EQUAL(size(Test4), 22);
   CHECK_EQUAL(norm_1(Test4), 124);
   CHECK_EQUAL(Test4[1], 8);
   CHECK_EQUAL(Test4[5], -6);

   Vector<double> Test5 = direct_sum(12U, range(2,6), 16U);
   CHECK_EQUAL(size(Test5), 6);
   CHECK_EQUAL(Test5[0], 12);
   CHECK_EQUAL(Test5[1], 2);
   CHECK_EQUAL(Test5[2], 3);
   CHECK_EQUAL(Test5[3], 4);
   CHECK_EQUAL(Test5[4], 5);
   CHECK_EQUAL(Test5[5], 16);

   Vector<int> X = range(20,30);
   Vector<double> Test6 = direct_sum(U, X);
   CHECK_EQUAL(size(Test6), 15);
   CHECK_EQUAL(norm_1(Test6), 265);

   Vector<std::string> S1(2);
   S1[0] = "hello";
   S1[1] = "world";
   Vector<std::string> S2(2);
   S2[0] = "goodbye";
   S2[1] = "world";

   Vector<std::string> S3 = direct_sum(S1,S2);
   CHECK_EQUAL(size(S3), 4);
   CHECK_EQUAL(S3[0], "hello");
   CHECK_EQUAL(S3[1], "world");
   CHECK_EQUAL(S3[2], "goodbye");
   CHECK_EQUAL(S3[3], "world");
}
