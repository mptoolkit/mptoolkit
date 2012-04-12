// -*- C++ -*- $Id$

#include "tensor_eigen.h"
#include "regularize.h"
#include "linearalgebra/eigen.h"
#include "linearalgebra/matrix_utility.h"

namespace Tensor
{

IrredTensor<LinearAlgebra::Matrix<double>, VectorBasis, VectorBasis>
DiagonalizeHermitian(IrredTensor<LinearAlgebra::Matrix<double>, VectorBasis, VectorBasis>& x)
{
   typedef IrredTensor<LinearAlgebra::Matrix<double>, VectorBasis, VectorBasis>
      TensorType;

   DEBUG_CHECK(is_scalar(x.TransformsAs()));
   DEBUG_CHECK_EQUAL(x.Basis1(), x.Basis2());
   DEBUG_CHECK(is_regular_basis(x.Basis1()))(x.Basis1());

   TensorType Result(x);
   for (iterator<TensorType>::type I = iterate(Result); I; ++I)
   {
      for (inner_iterator<TensorType>::type J = iterate(I); J; ++J)
      {
         DEBUG_CHECK_EQUAL(J.index1(), J.index2());
         LinearAlgebra::Vector<double> EVal = LinearAlgebra::DiagonalizeHermitian(*J);
         x(J.index1(), J.index2()) = LinearAlgebra::diagonal_matrix(EVal);
      }
   }

   return Result;
}


IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis>
DiagonalizeHermitian(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, 
                     VectorBasis, VectorBasis>& x)
{
   typedef IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis>
      TensorType;

   DEBUG_CHECK(is_scalar(x.TransformsAs()));
   DEBUG_CHECK_EQUAL(x.Basis1(), x.Basis2());
   DEBUG_CHECK(is_regular_basis(x.Basis1()))(x.Basis1());

   TensorType Result(x);
   for (iterator<TensorType>::type I = iterate(Result); I; ++I)
   {
      for (inner_iterator<TensorType>::type J = iterate(I); J; ++J)
      {
         DEBUG_CHECK_EQUAL(J.index1(), J.index2());
         LinearAlgebra::Vector<double> EVal = LinearAlgebra::DiagonalizeHermitian(*J);
         x(J.index1(), J.index2()) = LinearAlgebra::diagonal_matrix(EVal);
      }
   }

   return Result;
}

LinearAlgebra::Vector<double> 
EigenvaluesHermitian(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, 
                     VectorBasis, VectorBasis> const& x)
{
   LinearAlgebra::Vector<double> Result;
   for (std::size_t i = 0; i < x.Basis1().size(); ++i)
   {
      if (iterate_at(x.data(), i, i))
      {
         Result = direct_sum(Result, LinearAlgebra::EigenvaluesHermitian(x(i,i)));
      }
   }
   return Result;
}

void
InvertHPD(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis>& x)
{
   DEBUG_CHECK(is_scalar(x.TransformsAs()));
   DEBUG_CHECK_EQUAL(x.Basis1(), x.Basis2());
   DEBUG_CHECK(is_regular_basis(x.Basis1()))(x.Basis1());

   // we deliberately index the matrix here: if any components are missing
   // then we have a precondition failure (the matrix would be singular).
   for (std::size_t i = 0; i < x.Basis1().size(); ++i)
   {
      InvertHPD(x(i,i));
   }
}

void
InvertGeneral(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis>& x)
{
   DEBUG_CHECK(is_scalar(x.TransformsAs()));
   DEBUG_CHECK_EQUAL(x.Basis1(), x.Basis2());
   DEBUG_CHECK(is_regular_basis(x.Basis1()))(x.Basis1());

   // we deliberately index the matrix here: if any components are missing
   // then we have a precondition failure (the matrix would be singular).
   for (std::size_t i = 0; i < x.Basis1().size(); ++i)
   {
      InvertGeneral(x(i,i));
   }
}

void
InvertIrregularHPD(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis>& x)
{
   DEBUG_CHECK(is_scalar(x.TransformsAs()));
   DEBUG_CHECK_EQUAL(x.Basis1(), x.Basis2());
   IrredTensor<LinearAlgebra::Matrix<double>, VectorBasis, VectorBasis> U = Regularize(x.Basis1());
   x = triple_prod(U, x, herm(U));
   InvertHPD(x);
   x = triple_prod(herm(U), x, U);
}

void 
SingularValueDecompositionRegular(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, 
				  VectorBasis, VectorBasis> const& m, 
				  IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, 
				  VectorBasis, VectorBasis>& U,
                                  IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, 
					      VectorBasis, VectorBasis>& D,
				     IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, 
						 VectorBasis, VectorBasis>& Vh)
{
   typedef IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, 
      VectorBasis, VectorBasis> IrredT;
   DEBUG_CHECK(is_scalar(m.TransformsAs()));
   DEBUG_CHECK(is_regular_basis(m.Basis1()))(m.Basis1());
   DEBUG_CHECK(is_regular_basis(m.Basis2()))(m.Basis2());
   // make the basis for the D matrix
   std::vector<int> BasisMap(m.Basis1().size(), 0);
   VectorBasis DBasis(m.GetSymmetryList());
   for (unsigned i = 0; i < m.Basis1().size(); ++i)
   {
      for (unsigned j = 0; j < m.Basis2().size(); ++j)
      {
	 if (m.Basis1()[i] != m.Basis2()[j])
	    continue;

	 BasisMap[i] = DBasis.size();
	 DBasis.push_back(m.Basis1()[i], std::min(m.Basis1().dim(i), m.Basis2().dim(j)));
      }
   }

   U = IrredT(m.Basis1(), DBasis);
   D = IrredT(DBasis, DBasis);
   Vh = IrredT(DBasis, m.Basis2());

   for (unsigned i = 0; i < m.Basis1().size(); ++i)
   {
      for (unsigned j = 0; j < m.Basis2().size(); ++j)
      {
	 if (m.Basis1()[i] != m.Basis2()[j])
	    continue;

	 if (!iterate_at(m.data(), i, j))
	    continue;
	 set_element(U.data(), i,BasisMap[i], LinearAlgebra::Matrix<std::complex<double> >());
	 set_element(Vh.data(), BasisMap[i], j, LinearAlgebra::Matrix<std::complex<double> >());
	 LinearAlgebra::Vector<double> Dvec;
	 LinearAlgebra::SingularValueDecomposition(m(i,j), 
						   U(i,BasisMap[i]),
						   Dvec,
						   Vh(BasisMap[i], j));
	 set_element(D, BasisMap[i], BasisMap[i], diagonal_matrix(Dvec));
      }
   }
}

void 
SingularValueDecomposition(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, 
			               VectorBasis, VectorBasis> const& m, 
			   IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, 
			               VectorBasis, VectorBasis>& U,
			   IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, 
			               VectorBasis, VectorBasis>& D,
			   IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, 
			   VectorBasis, VectorBasis>& Vh)
{
   typedef IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, 
			               VectorBasis, VectorBasis> MatrixOperator;

   // If the basis is already regular, don't bother with the extra work
   if (is_regular_basis(m.Basis1()) && is_regular_basis(m.Basis2()))
   {
      SingularValueDecompositionRegular(m, U, D, Vh);
      return;
   }
   // else

   MatrixOperator U1 = Regularize(m.Basis1());
   MatrixOperator U2 = Regularize(m.Basis2());

   SingularValueDecompositionRegular(triple_prod(U1, m, herm(U2)), U, D, Vh);
   U = herm(U1) * U;
   Vh = Vh * U2;
}

LinearAlgebra::Vector<double>
SingularValues(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, 
	       VectorBasis, VectorBasis> const& m)
{
   typedef IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis> ITensor;

   ITensor U, D, Vh;
   SingularValueDecomposition(m, U, D, Vh);
   LinearAlgebra::Vector<double> Result(D.Basis1().total_dimension());
   int k = 0;
   for (unsigned i = 0; i < D.Basis1().size(); ++i)
   {
      for (int j = 0; j < D.Basis1().dim(i); ++j)
      {
	 if (iterate_at(D.data(), i, i))
	    Result[k++] = D(i,i)(j,j).real();
	 else
	    Result[k++] = 0.0;
      }
   }
   return Result;
}

void 
SingularValueDecomposition(IrredTensor<std::complex<double>, BasisList, BasisList> const& m, 
			   IrredTensor<std::complex<double>, BasisList, BasisList>& U,
			   IrredTensor<std::complex<double>, BasisList, BasisList>& D,
			   IrredTensor<std::complex<double>, BasisList, BasisList>& Vh)
{
   typedef IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, 
			               VectorBasis, VectorBasis> MatrixOperator;

   typedef IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, BasisList> MixedOperator;

   MixedOperator U1 = Regularize(m.Basis1());
   MixedOperator U2 = Regularize(m.Basis2());

   MatrixOperator UMatrix, DMatrix, VhMatrix;
   SingularValueDecompositionRegular(triple_prod(U1, m, herm(U2)), UMatrix, DMatrix, VhMatrix);

   MixedOperator Splitter = SplitBasis(DMatrix.Basis1());

   U = map_1x1_operator(triple_prod(herm(U1), UMatrix, Splitter));
   D = map_1x1_operator(triple_prod(herm(Splitter), DMatrix, Splitter));
   Vh = map_1x1_operator(triple_prod(herm(Splitter), VhMatrix, U2));
   return;
}

IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis>
InvertDiagonal(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, 
	       VectorBasis, VectorBasis> const& Op, double Cutoff)
{
   IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis>
      Result(Op.Basis2(), Op.Basis1());

   for (unsigned i = 0; i < Op.Basis1().size(); ++i)
   {
      //TRACE(i);
      if (iterate_at(Op.data(), i, i))
      {
	 Result(i,i) = LinearAlgebra::Matrix<double>(Op.Basis2().dim(i), Op.Basis1().dim(i), 0.0);
	 for (int j = 0; j < std::min(Op.Basis1().dim(i), Op.Basis2().dim(i)); ++j)
	 {
	    std::complex<double> x = Op(i,i)(j,j);
	    if (LinearAlgebra::norm_frob(x) < Cutoff) 
	       x = 0;
	    else 
	       x = 1.0 / x;
	    Result(i,i)(j,j) = x;
	 }
      }
      //TRACE(Result(i,i));
   }
   
   return Result;
}

IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis>
InvertDiagonal(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, 
	       VectorBasis, VectorBasis> const& Op)
{
   return InvertDiagonal(Op, std::sqrt(std::numeric_limits<double>::epsilon()));
}

IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis>
SqrtDiagonal(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, 
	       VectorBasis, VectorBasis> const& Op, double Tol)
{
   IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis>
      Result(Op.Basis2(), Op.Basis1());

   for (unsigned i = 0; i < Op.Basis1().size(); ++i)
   {
      if (iterate_at(Op.data(), i, i))
      {
	 Result(i,i) = LinearAlgebra::Matrix<double>(Op.Basis2().dim(i), Op.Basis1().dim(i), 0.0);
	 for (int j = 0; j < std::min(Op.Basis1().dim(i), Op.Basis2().dim(i)); ++j)
	 {
	    std::complex<double> x = Op(i,i)(j,j);
	    CHECK(LinearAlgebra::norm_frob(x.imag()) <= Tol);
	    if (x.real() <= 0)
	    {
	       CHECK(norm_frob(x.real()) <= Tol)(x.real())(Tol)(Op);
	       Result(i,i)(j,j) = 0.0;
	    }
	    else
	       Result(i,i)(j,j) = std::sqrt(x.real());
	 }
      }
      //TRACE(Result(i,i));
   }
   
   return Result;
}

// CholeskyFactorizeUpper

IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis>
CholeskyFactorizeUpperRegular(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, 
			      VectorBasis, VectorBasis> const& m)
{
   DEBUG_PRECONDITION_EQUAL(m.Basis1(), m.Basis2());

   typedef IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, 
      VectorBasis, VectorBasis> TensorType;

   TensorType Result = m;

   for (iterator<TensorType>::type I = iterate(Result); I; ++I)
   {
      for (inner_iterator<TensorType>::type J = iterate(I); J; ++J)
      {
	 LinearAlgebra::Matrix<std::complex<double> > x = *J;
	 CholeskyFactorizeUpper(x);
	 // zero out the lower triangular part
	 for (unsigned i = 1; i < size1(x); ++i)
	 {
	    for (unsigned j = 0; j < i; ++j)
	    {
	       x(i,j) = 0.0;
	    }
	 }
	 *J = x;
      }
   }
   return Result;
}

IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis>
CholeskyFactorizeUpper(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, 
		       VectorBasis, VectorBasis> const& m)
{
   DEBUG_PRECONDITION_EQUAL(m.Basis1(), m.Basis2());

   if (is_regular_basis(m.Basis1()))
      return CholeskyFactorizeUpperRegular(m);
   // else
   IrredTensor<LinearAlgebra::Matrix<double>, VectorBasis, VectorBasis> R = Regularize(m.Basis1());
   return triple_prod(herm(R), CholeskyFactorizeUpperRegular(triple_prod(R, m, herm(R))), R);
}

// CholeskyFactorizeLower

IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis>
CholeskyFactorizeLowerRegular(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, 
			      VectorBasis, VectorBasis> const& m)
{
   DEBUG_PRECONDITION_EQUAL(m.Basis1(), m.Basis2());

   typedef IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, 
      VectorBasis, VectorBasis> TensorType;

   TensorType Result = m;

   for (iterator<TensorType>::type I = iterate(Result); I; ++I)
   {
      for (inner_iterator<TensorType>::type J = iterate(I); J; ++J)
      {
	 LinearAlgebra::Matrix<std::complex<double> > x = *J;
	 CholeskyFactorizeLower(x);
	 // zero out the upper triangular part
	 unsigned const s2 = size2(x);
	 for (unsigned i = 0; i < size1(x); ++i)
	 {
	    for (unsigned j = i+1; j < s2; ++j)
	    {
	       x(i,j) = 0.0;
	    }
	 }
	 *J = x;
      }
   }
   return Result;
}

IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis>
CholeskyFactorizeLower(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, 
		       VectorBasis, VectorBasis> const& m)
{
   DEBUG_PRECONDITION_EQUAL(m.Basis1(), m.Basis2());

   if (is_regular_basis(m.Basis1()))
      return CholeskyFactorizeLowerRegular(m);
   // else
   IrredTensor<LinearAlgebra::Matrix<double>, VectorBasis, VectorBasis> R = Regularize(m.Basis1());
   return triple_prod(herm(R), CholeskyFactorizeLowerRegular(triple_prod(R, m, herm(R))), R);
}

IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis>
SingularFactorizeRegular(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, 
			 VectorBasis, VectorBasis> const& m)
{
   DEBUG_PRECONDITION_EQUAL(m.Basis1(), m.Basis2());

   typedef IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, 
      VectorBasis, VectorBasis> TensorType;

   TensorType Result = m;

   for (iterator<TensorType>::type I = iterate(Result); I; ++I)
   {
      for (inner_iterator<TensorType>::type J = iterate(I); J; ++J)
      {
	 *J = SingularFactorize(*J);
      }
   }
   return Result;
}

IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis>
SingularFactorize(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, 
		  VectorBasis, VectorBasis> const& m)
{
   DEBUG_PRECONDITION_EQUAL(m.Basis1(), m.Basis2());

   if (is_regular_basis(m.Basis1()))
      return SingularFactorizeRegular(m);
   // else
   IrredTensor<LinearAlgebra::Matrix<double>, VectorBasis, VectorBasis> R = Regularize(m.Basis1());
   return triple_prod(herm(R), SingularFactorize(triple_prod(R, m, herm(R))), R);
}

} // namespace Tensor
