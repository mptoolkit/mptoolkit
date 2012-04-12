// -*- C++ -*- $Id: tensor_eigen.cpp 809 2008-01-07 07:54:39Z ianmcc $

#include "tensor_exponential.h"
#include "linearalgebra/exponential.h"

namespace Tensor
{

IrredTensor<std::complex<double>, BasisList, BasisList>
Exponentiate(IrredTensor<std::complex<double>, BasisList, BasisList> const& m)
{
   typedef IrredTensor<std::complex<double>, BasisList, BasisList> TensorType;
   PRECONDITION(is_scalar(m.TransformsAs()))("Can only exponentiate a scalar operator")(m.TransformsAs());
   PRECONDITION_EQUAL(m.Basis1(), m.Basis2());

   using QuantumNumbers::QuantumNumber;

   TensorType Result(m.Basis1(), m.Basis2());

   // enumerate the quantum numbers in m
   std::set<QuantumNumber> UsedQ = QuantumNumbersInBasis(m.Basis1());

   // linearize the basis
   for (std::set<QuantumNumber>::const_iterator Q = UsedQ.begin(); Q != UsedQ.end(); ++Q)
   {
      std::map<int, int> LinearMap = LinearizeQuantumNumberSubspace(m.Basis1(), *Q);
      LinearAlgebra::Matrix<std::complex<double> > M(LinearMap.size(), LinearMap.size(), 0.0);
      for (const_iterator<TensorType>::type I = iterate(m); I; ++I)
      {
	 if (m.Basis1()[I.index()] != *Q)
	    continue;
	 for (const_inner_iterator<TensorType>::type J = iterate(I); J; ++J)
	 {
	    if (m.Basis2()[J.index2()] != *Q)
	       continue;

	    M(LinearMap[J.index1()], LinearMap[J.index2()]) = *J;
	 }
      }

      M = LinearAlgebra::Exponentiate(1.0, M);
      for (std::map<int, int>::const_iterator I = LinearMap.begin(); I != LinearMap.end(); ++I)
      {
	 for (std::map<int, int>::const_iterator J = LinearMap.begin(); J != LinearMap.end(); ++J)
	 {
	    std::complex<double> x = M(I->second, J->second);
	    if (LinearAlgebra::norm_frob(x) > 1e-14)
	       Result(I->first, J->first) = x;
	 }
      }
   }

   return Result;
}

} // namespace Tensor
