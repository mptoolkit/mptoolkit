// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/deductiontest.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
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
