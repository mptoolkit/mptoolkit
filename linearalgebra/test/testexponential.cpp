// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/test/testexponential.cpp
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
//
// test the dense matrix exponential function.
// The main thing we want to do is see if it works for small inputs.

#include "linearalgebra/exponential.h"

#include <cmath>
#include <complex>
#include <limits>

using namespace LinearAlgebra;

namespace
{

using Complex = std::complex<double>;

void CheckClose(Complex Actual, Complex Expected, double Tolerance)
{
   CHECK_COMPARE(std::abs(Actual - Expected), <, Tolerance)(Actual)(Expected);
}

void CheckFiniteBounded(Complex Value, double Bound)
{
   CHECK(std::isfinite(Value.real()))(Value);
   CHECK(std::isfinite(Value.imag()))(Value);
   CHECK_COMPARE(std::abs(Value), <, Bound)(Value);
}

void CheckRelativeClose(Complex Actual, Complex Expected, double Tolerance)
{
   CHECK(std::isfinite(Actual.real()))(Actual);
   CHECK(std::isfinite(Actual.imag()))(Actual);
   CHECK_COMPARE(std::abs((Actual / Expected) - Complex(1.0)), <, Tolerance)(Actual)(Expected);
}

} // namespace

int main()
{
   double const Tolerance = 1e-10;

   {
      Matrix<Complex> M(3, 3, 0.0);
      M(0,0) = 0.25;
      M(1,1) = -0.5;
      M(2,2) = 1.0;

      Matrix<Complex> Result = Exponentiate(2.0, M);
      CheckClose(Result(0,0), std::exp(0.5), Tolerance);
      CheckClose(Result(1,1), std::exp(-1.0), Tolerance);
      CheckClose(Result(2,2), std::exp(2.0), Tolerance);
      CheckClose(Result(0,1), 0.0, Tolerance);
   }

   {
      Matrix<Complex> M(2, 2, 0.0);
      M(0,1) = 2.0;

      Matrix<Complex> Result = Exponentiate(1.0, M);
      CheckClose(Result(0,0), 1.0, Tolerance);
      CheckClose(Result(0,1), 2.0, Tolerance);
      CheckClose(Result(1,0), 0.0, Tolerance);
      CheckClose(Result(1,1), 1.0, Tolerance);
   }

   {
      double const Theta = 0.25;
      Matrix<Complex> M(2, 2, 0.0);
      M(0,1) = -Theta;
      M(1,0) = Theta;

      Matrix<Complex> Result = Exponentiate(1.0, M);
      CheckClose(Result(0,0), std::cos(Theta), Tolerance);
      CheckClose(Result(0,1), -std::sin(Theta), Tolerance);
      CheckClose(Result(1,0), std::sin(Theta), Tolerance);
      CheckClose(Result(1,1), std::cos(Theta), Tolerance);
   }

   {
      double const Theta = 1e80;
      Matrix<Complex> M(2, 2, 0.0);
      M(0,1) = -Theta;
      M(1,0) = Theta;

      Matrix<Complex> Result = Exponentiate(1.0, M);
      CheckFiniteBounded(Result(0,0), 2.0);
      CheckFiniteBounded(Result(0,1), 2.0);
      CheckFiniteBounded(Result(1,0), 2.0);
      CheckFiniteBounded(Result(1,1), 2.0);
   }

   {
      double const Value = 0.75 * std::numeric_limits<double>::max();
      Matrix<Complex> M(3, 3, 0.0);
      M(0,2) = Value;
      M(1,2) = -Value;

      Matrix<Complex> Result = Exponentiate(1.0, M);
      CheckClose(Result(0,0), 1.0, Tolerance);
      CheckClose(Result(1,1), 1.0, Tolerance);
      CheckClose(Result(2,2), 1.0, Tolerance);
      CheckClose(Result(0,1), 0.0, Tolerance);
      CheckClose(Result(1,0), 0.0, Tolerance);
      CheckRelativeClose(Result(0,2), Value, Tolerance);
      CheckRelativeClose(Result(1,2), -Value, Tolerance);
   }
}
