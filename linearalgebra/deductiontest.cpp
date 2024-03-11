// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/deductiontest.cpp
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

#include <iostream>
#include <iomanip>

using namespace std;

template <typename T>
struct base {};

template <typename T>
struct derived : public base<T> {};

typedef char yes[2];
typedef char no[1];

no& foo(...);

template <typename T>
yes& foo(base<T> const& b);

int main()
{
  int i;
  base<long> b;
  derived<double> d;
  cout << sizeof(foo(i)) << '\n' << sizeof(foo(b)) << '\n' << sizeof(foo(d)) << '\n';
}
