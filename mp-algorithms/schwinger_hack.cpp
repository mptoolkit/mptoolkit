// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/schwinger_hack.cpp
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

#include "schwinger_hack.h"
#include "common/environment.h"

// Schwinger model hacks
bool HackSchwinger_E = true;
bool HackSchwinger_F = true;
double HackSchwinger_Field = getenv_or_default("MPE", 1.0);
bool GaugeFlip = false;
//bool GaugeFlip = true;

MatrixOperator
GetQ(MatrixOperator LeftIdentity, MatrixOperator const& RightIdentity, int Power)
{
   VectorBasis b = LeftIdentity.Basis1();  // Basis1() and Basis2() are the same here
   MatrixOperator Q(b, b);
   for (int i = 0; i < b.size(); ++i)
   {
      Q(i,i) = std::pow(casimir(b[i],0), Power) * LeftIdentity(i,i);
   }
   return Q;
}

// Returns the expectation value of the first casimir operator of the basis,
// optionally at some power
std::complex<double>
GetQuantumNumberExpectation(MatrixOperator LeftIdentity, MatrixOperator const& RightIdentity, int Power)
{
   return inner_prod(GetQ(LeftIdentity, RightIdentity, Power), RightIdentity);
}
