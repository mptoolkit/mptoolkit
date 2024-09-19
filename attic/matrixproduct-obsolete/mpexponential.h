// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// attic/matrixproduct-obsolete/mpexponential.h
//
// Copyright (C) 2004-2017 Ian McCulloch <ian@qusim.net>
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

#if !defined(MPEXPONENTIAL_H_JHDIO48UY849Y894YWP89EYH)
#define MPEXPONENTIAL_H_JHDIO48UY849Y894YWP89EYH

#include "splitoperator.h"

std::pair<MPOpComponent, MPOpComponent>
TwoSiteExponential(SimpleOperator const& A, SimpleOperator const& B, std::complex<double> x);

// calculates exp(x * (A \otimes B))
// where (A \otimes B) is a scalar operator.
void TwoSiteExponential(MPOpComponent& A, MPOpComponent& B, std::complex<double> x);

// calculates the exponential of a bond operator.  Op must have all except
// two adjacent terms proportional to the identity.  Op must transform as a scalar.
SplitOperator BondExponential(std::complex<double> x, SplitOperator const& Op);

#endif
