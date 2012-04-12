// -*- C++ -*- $Id$
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
