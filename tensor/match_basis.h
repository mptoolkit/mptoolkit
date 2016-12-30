// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// tensor/match_basis.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
//
// Constructs a bi-projection to match two
// basis sets (eg, for PWFRG)

#if !defined(MATCH_BASIS_DSH458743Y67YE)
#define MATCH_BASIS_DSH458743Y67YE

#include "basis.h"
#include "tensor.h"

namespace Tensor
{

IrredTensor<LinearAlgebra::Matrix<double>, VectorBasis, VectorBasis>
MatchBasis(VectorBasis const& B1, VectorBasis const& B2);

   // match basis starting from the last entry and working backwards
IrredTensor<LinearAlgebra::Matrix<double>, VectorBasis, VectorBasis>
MatchBasisReverse(VectorBasis const& B1, VectorBasis const& B2);

} // namespace Tensor

#endif
