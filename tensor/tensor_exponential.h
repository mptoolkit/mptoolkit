// -*- C++ -*- $Id: tensor_eigen.h 809 2008-01-07 07:54:39Z ianmcc $

#if !defined(TENSOR_EXPONENTIAL_H_DHCKJDHUREYT7845Y78Y78TY78TY78T)
#define TENSOR_EXPONENTIAL_H_DHCKJDHUREYT7845Y78Y78TY78TY78T

#include "tensor.h"

namespace Tensor
{

// exponentiate a scalar operator
IrredTensor<std::complex<double>, BasisList, BasisList>
Exponentiate(IrredTensor<std::complex<double>, BasisList, BasisList> const& m);

} // namespace Tensor

#endif
