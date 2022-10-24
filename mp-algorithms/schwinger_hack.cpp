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
bool HackSchwinger_E = false;
bool HackSchwinger_F = false;
bool GaugeFlip = false;
//bool GaugeFlip = true;

double HackSchwinger_Field = 0;
double HackSchwinger_Bg = 0;

void SetSchwingerFieldFromAttributes(AttributeList const& A)
{
   if (!A["SchwingerField"].empty())
   {
      HackSchwinger_E = true;
      HackSchwinger_F = true;
      HackSchwinger_Field = A["SchwingerField"].as<double>();
      std::cerr << "Schwinger field is " << HackSchwinger_Field << '\n';
   }
   if (!A["SchwingerBackground"].empty())
   {
      HackSchwinger_Bg = A["SchwingerBackground"].as<double>();
      std::cerr << "Schwinger background field is " << HackSchwinger_Bg << '\n';
   }
}

MatrixOperator
GetSchwingerBoundaryMatrix(MatrixOperator LeftIdentity, MatrixOperator const& RightIdentity, int Power)
{
   VectorBasis b = LeftIdentity.Basis1();  // Basis1() and Basis2() are the same here
   MatrixOperator Q(b, b);
   for (int i = 0; i < b.size(); ++i)
   {
      Q(i,i) = std::pow(casimir(b[i],0) + HackSchwinger_Bg, Power) * LeftIdentity(i,i);
   }
   return Q;
}

// Returns the expectation value of the first casimir operator of the basis,
// optionally at some power
std::complex<double>
GetSchwingerBoundaryCorrection(MatrixOperator LeftIdentity, MatrixOperator const& RightIdentity, int Power)
{
   return inner_prod(GetSchwingerBoundaryMatrix(LeftIdentity, RightIdentity, Power), RightIdentity);
}
