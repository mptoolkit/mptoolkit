// -*- C++ -*- $Id$

namespace Tensor
{

template <typename T>
IrredTensor<LinearAlgebra::DiagonalMatrix<T>, VectorBasis, VectorBasis>
InvertDiagonal(IrredTensor<LinearAlgebra::DiagonalMatrix<T>, 
	       VectorBasis, VectorBasis> const& Op, double Cutoff)
{
   IrredTensor<LinearAlgebra::DiagonalMatrix<T>, VectorBasis, VectorBasis>
      Result(Op.Basis2(), Op.Basis1());

   for (unsigned i = 0; i < Op.Basis1().size(); ++i)
   {
      if (iterate_at(Op.data(), i, i))
      {
	 Result(i,i) = LinearAlgebra::DiagonalMatrix<T>(Op.Basis2().dim(i), Op.Basis1().dim(i), 0.0);
	 for (int j = 0; j < std::min(Op.Basis1().dim(i), Op.Basis2().dim(i)); ++j)
	 {
	    T x = Op(i,i)(j,j);
	    if (LinearAlgebra::norm_frob(x) < Cutoff) 
	       x = 0;
	    else 
	       x = 1.0 / x;
	    Result(i,i)(j,j) = x;
	 }
      }
   }
   return Result;
}

template <typename T>
IrredTensor<LinearAlgebra::DiagonalMatrix<T>, VectorBasis, VectorBasis>
InvertDiagonal(IrredTensor<LinearAlgebra::DiagonalMatrix<T>, 
	       VectorBasis, VectorBasis> const& Op)
{
   return InvertDiagonal(Op, std::sqrt(std::numeric_limits<double>::epsilon()));
}

} // namespace Tensor
