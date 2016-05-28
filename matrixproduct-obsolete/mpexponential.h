// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// matrixproduct-obsolete/mpexponential.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
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
