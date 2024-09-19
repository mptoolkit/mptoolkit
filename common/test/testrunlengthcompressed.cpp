// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// common/test/testrunlengthcompressed.cpp
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

#include "common/runlengthcompressed.h"
#include <iostream>

template <typename T>
void print(run_length_compressed<T> const& x)
{
   std::copy(x.begin(), x.end(), std::ostream_iterator<T>(std::cout, "\n"));
}

template <typename T>
void print_reverse(run_length_compressed<T> const& x)
{
   std::reverse_copy(x.begin(), x.end(), std::ostream_iterator<T>(std::cout, "\n"));
}

int main()
{
   run_length_compressed<int> x;

   CHECK(x.empty());
   CHECK_EQUAL(x.size(), 0);

   std::vector<int> xTest(x.begin(), x.end());
   CHECK(xTest.empty());

   x.push_back(5);
   CHECK(!x.empty());
   CHECK_EQUAL(x.size(), 1);

   xTest.assign(x.begin(), x.end());
   CHECK_EQUAL(xTest.size(), 1);
   CHECK_EQUAL(xTest[0], 5);

   x = repeat(x, 3);
   CHECK_EQUAL(x.size(), 3);
   xTest.assign(x.begin(), x.end());
   CHECK_EQUAL(xTest.size(), 3);
   CHECK_EQUAL(xTest[0], 5);
   CHECK_EQUAL(xTest[1], 5);
   CHECK_EQUAL(xTest[2], 5);

   run_length_compressed<int> y;
   y.push_back(7);
   x.push_back(y);
   CHECK_EQUAL(x.size(), 4);

   x = repeat(x, 4);
   CHECK_EQUAL(x.node_count(), 3);  // 3 nodes: the repeat of 3x'5', the array of r=7,555 and the repeat of 4xr.
   CHECK_EQUAL(x.leaf_count(), 2);  // 2 leaves: 5 and 7
   CHECK_EQUAL(x.size(), 16);
   CHECK_EQUAL(x.front(), 5);
   CHECK_EQUAL(x.back(), 7);

   // make sure that adding a empty node does nothing
   x.push_back(run_length_repeat<int>(0, run_length_compressed<int>(4)));
   CHECK_EQUAL(x.size(), 16);
   CHECK_EQUAL(x.node_count(), 3);
   CHECK_EQUAL(x.leaf_count(), 2);

   x.push_back(run_length_array<int>());
   CHECK_EQUAL(x.size(), 16);
   CHECK_EQUAL(x.node_count(), 3);
   CHECK_EQUAL(x.leaf_count(), 2);

   x.push_front(run_length_repeat<int>(0, run_length_compressed<int>(4)));
   CHECK_EQUAL(x.size(), 16);
   CHECK_EQUAL(x.node_count(), 3);
   CHECK_EQUAL(x.leaf_count(), 2);

   x.push_front(run_length_array<int>());
   CHECK_EQUAL(x.size(), 16);
   CHECK_EQUAL(x.node_count(), 3);
   CHECK_EQUAL(x.leaf_count(), 2);

   xTest.clear();
   std::copy(x.begin(), x.end(), std::back_insert_iterator<std::vector<int> >(xTest));
   CHECK_EQUAL(xTest.size(), 16);
   int xCheck[] = {5,5,5,7,5,5,5,7,5,5,5,7,5,5,5,7};
   CHECK(xTest == std::vector<int>(xCheck,xCheck+16));

   xTest.clear();
   std::reverse_copy(x.begin(), x.end(), std::back_insert_iterator<std::vector<int> >(xTest));
   CHECK_EQUAL(xTest.size(), 16);
   int xCheckRev[] = {7,5,5,5,7,5,5,5,7,5,5,5,7,5,5,5};
   CHECK(xTest == std::vector<int>(xCheckRev,xCheckRev+16));

   run_length_compressed<int>::const_iterator I = x.begin();
   std::advance(I, 6);
   xTest.assign(x.begin(), I);
   CHECK_EQUAL(xTest.size(), 6);
   CHECK(xTest == std::vector<int>(xCheck,xCheck+xTest.size()));

   xTest.assign(I, x.end());
   CHECK_EQUAL(xTest.size(), 10);
   CHECK(xTest == std::vector<int>(xCheck+6,xCheck+16));

   std::advance(I, 2);
   xTest.assign(x.begin(), I);
   CHECK_EQUAL(xTest.size(), 8);
   CHECK(xTest == std::vector<int>(xCheck,xCheck+xTest.size()));

   xTest.assign(I, x.end());
   CHECK_EQUAL(xTest.size(), 8);
   CHECK(xTest == std::vector<int>(xCheck+8,xCheck+16));
}
