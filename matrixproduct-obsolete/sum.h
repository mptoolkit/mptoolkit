// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// matrixproduct-obsolete/sum.h
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

#include "matrixproduct.h"

template <typename Component>
MatrixProduct<Component>
fancy_sum(MatrixProduct<Component> A, MatrixProduct<Component> B, int NumStates)
{
   typedef typename Component::OperatorType OperatorType;

   if (A.is_null()) return B;
   if (B.is_null()) return A;

   QuantumNumber Ident(A.GetSymmetryList());  // the scalar quantum number

   // We start from the right hand side here...
   while (A.RightSize() > 1) A.RotateRight();
   while (B.RightSize() > 1) B.RotateRight();
   
