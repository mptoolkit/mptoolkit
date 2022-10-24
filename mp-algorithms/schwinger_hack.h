// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/schwinger_hack.h
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

#if !defined(MPTOOLKIT_MP_ALGORITHMS_SCHWINGER_HACK_H)
#define MPTOOLKIT_MP_ALGORITHMS_SCHWINGER_HACK_H

#include "triangular_mpo_solver.h"
#include "interface/attributes.h"

extern bool HackSchwinger_E;
extern bool HackSchwinger_F;
extern double HackSchwinger_Field;
extern double HackSchwinger_Bg;
extern bool GaugeFlip;

void SetSchwingerFieldFromAttributes(AttributeList const& A);

MatrixOperator SchwingerBoundaryMatrix(MatrixOperator LeftIdentity, MatrixOperator const& RightIdentity, int Power = 1);

// Returns the expectation value of the first casimir operator of the basis,
// optionally at some power
std::complex<double>
GetSchwingerBoundaryCorrection(MatrixOperator LeftIdentity, MatrixOperator const& RightIdentity, int Power = 1);

#endif
