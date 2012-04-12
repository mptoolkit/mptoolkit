// -*- C++ -*- $Id$

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
