// -*- C++ -*- $Id$

namespace Tensor
{

template <typename T>
IrredTensor<T, BasisList, BasisList>
map_1x1_operator(IrredTensor<LinearAlgebra::Matrix<T>, BasisList, BasisList> const& Op)
{
   typedef IrredTensor<LinearAlgebra::Matrix<T>, BasisList, BasisList> MatOpType;
   IrredTensor<T, BasisList, BasisList> Result(Op.Basis1(), Op.Basis2(), Op.TransformsAs());

   for (typename const_iterator<MatOpType>::type I = iterate(Op); I; ++I)
   {
      for (typename const_inner_iterator<MatOpType>::type J = iterate(I); J; ++J)
      {
         DEBUG_PRECONDITION_EQUAL(size1(*J), 1);
         DEBUG_PRECONDITION_EQUAL(size2(*J), 1);
         Result(J.index1(), J.index2()) = (*J)(0,0);
      }
   }
   return Result;
}

template <typename T>
IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis>
ToMatrixOperator(IrredTensor<T, BasisList, BasisList> const& Op)
{
   IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis> 
      Result(VectorBasis(Op.Basis1()), VectorBasis(Op.Basis2()), Op.TransformsAs());
   Result.data() = LinearAlgebra::scalar(LinearAlgebra::Matrix<T>(1,1,1.0)) * Op.data();
   return Result;
}

} // namespace Tensor
