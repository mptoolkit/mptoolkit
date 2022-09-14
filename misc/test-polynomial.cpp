// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// misc/test-polynomial.h
//
// Copyright (C) 2022 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "common/polynomial.h"
#include <complex>
#include <iostream>

using Poly = Polynomial<std::complex<double>>;


int main()
{
   Poly p;
   p[1] = 0.5;
   p[2] = 1.0;
   p[4] = 4.0;
   std::cout << p(1.0) << '\n';
   std::cout << p.derivative(1.0, 2) << '\n';
   // p(1) should be 0.5 + 1 + 4 = 5.5
   // derivative of p is 0.5 + 2x + 16x^3
   // second derivative of p is 2 + 48*x^2
   // second derivative evaluated at 1.0 is 50
}
