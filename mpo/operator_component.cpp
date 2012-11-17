// -*- C++ -*- $Id$

#include "operator_component.h"
#include "tensor/tensorproduct.h"
#include "linearalgebra/eigen.h"
#include "tensor/regularize.h"
#include "linearalgebra/matrix_utility.h"

//using namespace LinearAlgebra;
using LinearAlgebra::operator*;

OperatorComponent::OperatorComponent(BasisList const& LocalB, 
                                     BasisList const& B1, BasisList const& B2)
   : LocalBasis1_(LocalB), LocalBasis2_(LocalB), Basis1_(B1), Basis2_(B2),
     Data_(Basis1_.size(), Basis2_.size())
{
}

OperatorComponent::OperatorComponent(BasisList const& LocalB1, BasisList const& LocalB2, 
                                     BasisList const& B1, BasisList const& B2)
   : LocalBasis1_(LocalB1), LocalBasis2_(LocalB2), Basis1_(B1), Basis2_(B2),
     Data_(Basis1_.size(), Basis2_.size())
{
}

std::ostream&
operator<<(std::ostream& out, OperatorComponent const& op)
{
   out << "OperatorComponent: Basis1=" << op.Basis1() << "\nBasis2=" << op.Basis2() << '\n';
   for (LinearAlgebra::const_iterator<OperatorComponent>::type I = iterate(op); I; ++I)
   {
      for (LinearAlgebra::const_inner_iterator<OperatorComponent>::type J = iterate(I); J; ++J)
      {
         out << "element (" << J.index1() << "," << J.index2() << ") = " << *J << '\n';
      }
   }
   return out;
}

PStream::opstream& 
operator<<(PStream::opstream& out, OperatorComponent const& Op)
{
   out << Op.LocalBasis1_ << Op.LocalBasis2_ 
       << Op.Basis1_ << Op.Basis2_
       << Op.Data_
      ;
   return out;
}

PStream::ipstream& 
operator>>(PStream::ipstream& in, OperatorComponent& Op)
{
   in >> Op.LocalBasis1_ >> Op.LocalBasis2_ 
      >> Op.Basis1_ >> Op.Basis2_
      >> Op.Data_
      ;
   return in;
}

OperatorComponent 
operator+(OperatorComponent const& A, OperatorComponent const& B)
{
   DEBUG_CHECK_EQUAL(A.LocalBasis1(), B.LocalBasis1());
   DEBUG_CHECK_EQUAL(A.LocalBasis2(), B.LocalBasis2());
   DEBUG_CHECK_EQUAL(A.Basis1(), B.Basis1());
   DEBUG_CHECK_EQUAL(A.Basis2(), B.Basis2());
   OperatorComponent Result(A.LocalBasis1(), A.LocalBasis2(), A.Basis1(), A.Basis2());
   Result.data() = A.data() + B.data();
   return Result;
}

OperatorComponent
tensor_sum(OperatorComponent const& A, OperatorComponent const& B,
           SumBasis<BasisList> const& B1, SumBasis<BasisList> const& B2)
{
   DEBUG_CHECK_EQUAL(A.LocalBasis1(), B.LocalBasis1());
   DEBUG_CHECK_EQUAL(A.LocalBasis2(), B.LocalBasis2());

   OperatorComponent Result(A.LocalBasis1(), A.LocalBasis2(), B1.Basis(), B2.Basis());
   
   // set the matrix elements corresponding to operator A (subspace 0 in the SumBasis)
   for (LinearAlgebra::const_iterator<OperatorComponent>::type I = iterate(A); I; ++I)
   {
      for (LinearAlgebra::const_inner_iterator<OperatorComponent>::type J = iterate(I); J; ++J)
      {
         Result(B1(0,J.index1()), B2(0,J.index2())) = *J;
      }
   }

   // set the matrix elements corresponding to operator A (subspace 1 in the SumBasis)
   for (LinearAlgebra::const_iterator<OperatorComponent>::type I = iterate(B); I; ++I)
   {
      for (LinearAlgebra::const_inner_iterator<OperatorComponent>::type J = iterate(I); J; ++J)
      {
         Result(B1(1,J.index1()), B2(1,J.index2())) = *J;
      }
   }
   return Result;
}

OperatorComponent 
tensor_row_sum(OperatorComponent const& A, 
	       OperatorComponent const& B, 
	       SumBasis<BasisList> const& B2)
{
   DEBUG_CHECK_EQUAL(A.LocalBasis1(), B.LocalBasis1());
   DEBUG_CHECK_EQUAL(A.LocalBasis2(), B.LocalBasis2());
   DEBUG_CHECK_EQUAL(A.Basis1(), B.Basis1());

   OperatorComponent Result(A.LocalBasis1(), A.LocalBasis2(), A.Basis1(), B2.Basis());
   
   // set the matrix elements corresponding to operator A (subspace 0 in the SumBasis)
   for (LinearAlgebra::const_iterator<OperatorComponent>::type I = iterate(A); I; ++I)
   {
      for (LinearAlgebra::const_inner_iterator<OperatorComponent>::type J = iterate(I); J; ++J)
      {
         Result(J.index1(),B2(0, J.index2())) = *J;
      }
   }

   // set the matrix elements corresponding to operator A (subspace 1 in the SumBasis)
   for (LinearAlgebra::const_iterator<OperatorComponent>::type I = iterate(B); I; ++I)
   {
      for (LinearAlgebra::const_inner_iterator<OperatorComponent>::type J = iterate(I); J; ++J)
      {
         Result(J.index1(),B2(1, J.index2())) = *J;
      }
   }
   return Result;
}
 
OperatorComponent 
tensor_col_sum(OperatorComponent const& A, 
	       OperatorComponent const& B, 
	       SumBasis<BasisList> const& B1)
{
   DEBUG_CHECK_EQUAL(A.LocalBasis1(), B.LocalBasis1());
   DEBUG_CHECK_EQUAL(A.LocalBasis2(), B.LocalBasis2());
   DEBUG_CHECK_EQUAL(A.Basis2(), B.Basis2());

   OperatorComponent Result(A.LocalBasis1(), A.LocalBasis2(), B1.Basis(), A.Basis2());
   
   // set the matrix elements corresponding to operator A (subspace 0 in the SumBasis)
   for (LinearAlgebra::const_iterator<OperatorComponent>::type I = iterate(A); I; ++I)
   {
      for (LinearAlgebra::const_inner_iterator<OperatorComponent>::type J = iterate(I); J; ++J)
      {
         Result(B1(0,J.index1()), J.index2()) = *J;
      }
   }

   // set the matrix elements corresponding to operator A (subspace 1 in the SumBasis)
   for (LinearAlgebra::const_iterator<OperatorComponent>::type I = iterate(B); I; ++I)
   {
      for (LinearAlgebra::const_inner_iterator<OperatorComponent>::type J = iterate(I); J; ++J)
      {
         Result(B1(1,J.index1()), J.index2()) = *J;
      }
   }
   return Result;
}

OperatorComponent prod(OperatorComponent const& A, SimpleOperator const& Op)
{
   PRECONDITION(is_scalar(Op.TransformsAs()))(Op.TransformsAs());
   DEBUG_PRECONDITION_EQUAL(A.Basis2(), Op.Basis1());

   OperatorComponent Result(A.LocalBasis1(), A.LocalBasis2(), A.Basis1(), Op.Basis2());
   // for our convention on the coupling constants, this should be a trivial operation
   // of a standard matrix-matrix multiply (the nested operation is a SimpleRedOperator * scalar)
   Result.data() = A.data() * Op.data(); 
   return Result;
}

OperatorComponent prod(SimpleOperator const& Op, OperatorComponent const& A)
{
   PRECONDITION(is_scalar(Op.TransformsAs()))(Op.TransformsAs());
   DEBUG_PRECONDITION_EQUAL(Op.Basis2(), A.Basis1());

   OperatorComponent Result(A.LocalBasis1(), A.LocalBasis2(), Op.Basis1(), A.Basis2());
   Result.data() = Op.data() * A.data(); 
   return Result;
}

OperatorComponent prod(OperatorComponent const& A, HermitianProxy<SimpleOperator> const& Op)
{
   PRECONDITION(is_scalar(Op.base().TransformsAs()))(Op.base().TransformsAs());
   DEBUG_PRECONDITION_EQUAL(A.Basis2(), Op.base().Basis2());

   OperatorComponent Result(A.LocalBasis1(), A.LocalBasis2(), A.Basis1(), Op.base().Basis1());
   // for our convention on the coupling constants, this should be a trivial operation
   // of a standard matrix-matrix multiply (the nested operation is a SimpleRedOperator * scalar)
   Result.data() = A.data() * herm(Op.base().data()); 
   return Result;
}

OperatorComponent prod(HermitianProxy<SimpleOperator> const& Op, OperatorComponent const& A)
{
   PRECONDITION(is_scalar(Op.base().TransformsAs()))(Op.base().TransformsAs());
   DEBUG_PRECONDITION_EQUAL(Op.base().Basis1(), A.Basis1());

   OperatorComponent Result(A.LocalBasis1(), A.LocalBasis2(), Op.base().Basis2(), A.Basis2());
   Result.data() = herm(Op.base().data()) * A.data(); 
   return Result;
}

OperatorComponent 
triple_prod(SimpleOperator const& x, 
            OperatorComponent const& Op, 
            LinearAlgebra::HermitianProxy<SimpleOperator> const& y)
{
   return prod(x, prod(Op, y));
}

OperatorComponent local_tensor_prod(OperatorComponent const& A, OperatorComponent const& B)
{
   DEBUG_PRECONDITION_EQUAL(A.Basis2(), B.Basis1());

   ProductBasis<BasisList, BasisList> LB1(A.LocalBasis1(), B.LocalBasis1());
   ProductBasis<BasisList, BasisList> LB2(A.LocalBasis2(), B.LocalBasis2());

   OperatorComponent Result(LB1.Basis(), LB2.Basis(), A.Basis1(), B.Basis2());

   typedef OperatorComponent::data_type MatType;

   Result.data() = LinearAlgebra::MatrixMatrixMultiplication<MatType, MatType, Tensor::TensorProd<SimpleRedOperator, SimpleRedOperator> >()(A.data(), B.data(), Tensor::TensorProd<SimpleRedOperator, SimpleRedOperator>());

   return Result;
}

#if 0
//  *** This isn't correct for SU(2), the direct_product of the matrix doesn't
// give the same thing as the ProductBasis. ***
OperatorComponent aux_tensor_prod(OperatorComponent const& A, OperatorComponent const& B)
{
   DEBUG_PRECONDITION_EQUAL(A.LocalBasis2(), B.LocalBasis1());
   ProductBasis<BasisList, BasisList> B1(A.Basis1(), B.Basis1());
   ProductBasis<BasisList, BasisList> B2(A.Basis2(), B.Basis2());

   OperatorComponent Result(A.LocalBasis1(), B.LocalBasis2(), B1.Basis(), B2.Basis());
   Result.data() = direct_product(A.data(), B.data(), LinearAlgebra::Multiplication<SimpleRedOperator, SimpleRedOperator>());
   return Result;
}

//  *** This isn't correct for SU(2), the direct_product of the matrix doesn't
// give the same thing as the ProductBasis. ***
OperatorComponent global_tensor_prod(OperatorComponent const& A, OperatorComponent const& B)
{
   ProductBasis<BasisList, BasisList> B1(A.Basis1(), B.Basis1());
   ProductBasis<BasisList, BasisList> B2(A.Basis2(), B.Basis2());
   ProductBasis<BasisList, BasisList> LB1(A.LocalBasis1(), B.LocalBasis1());
   ProductBasis<BasisList, BasisList> LB2(A.LocalBasis2(), B.LocalBasis2());

   OperatorComponent Result(LB1.Basis(), LB2.Basis(), B1.Basis(), B2.Basis());
   Result.data() = direct_product(A.data(), B.data(), Tensor::TensorProd<SimpleRedOperator, SimpleRedOperator>());
   return Result;
}
#endif

OperatorComponent exchange(OperatorComponent const& A)
{
   // Swap the local and aux indices.  A(i,j)[k](r,s) -> Result(r,s)[k](i,j)
   OperatorComponent Result(A.Basis1(), A.Basis2(), A.LocalBasis1(), A.LocalBasis2());
   // Iterate over the components in A, first index
   for (LinearAlgebra::const_iterator<OperatorComponent>::type I = iterate(A); I; ++I)
   {
      // second index in A
      for (LinearAlgebra::const_inner_iterator<OperatorComponent>::type J = iterate(I); J; ++J)
      {
         // Iterate over the irreducible components of A(I,J)
         for (SimpleRedOperator::const_iterator k = J->begin(); k != J->end(); ++k)
         {
            // *k is an irreducible operator.  Iterate over the components of this operator
            for (LinearAlgebra::const_iterator<SimpleOperator>::type R = iterate(*k); R; ++R)
            {
               for (LinearAlgebra::const_inner_iterator<SimpleOperator>::type 
                       S = iterate(R); S; ++S)
               {
                  Result(S.index1(), S.index2()).project(k->TransformsAs())(J.index1(), J.index2()) = *S;
               }
            }
         }
      }
   }
   return Result;
}

StateComponent
operator_prod(OperatorComponent const& M,
              StateComponent const& A, 
              StateComponent const& E,
              HermitianProxy<StateComponent> const& B)
{
   DEBUG_PRECONDITION_EQUAL(M.LocalBasis1(), A.LocalBasis());
   DEBUG_PRECONDITION_EQUAL(M.LocalBasis2(), B.base().LocalBasis());
   DEBUG_PRECONDITION_EQUAL(M.Basis2(), E.LocalBasis());
   DEBUG_PRECONDITION_EQUAL(A.Basis2(), E.Basis1());
   DEBUG_PRECONDITION_EQUAL(E.Basis2(), B.base().Basis2());

   StateComponent Result(M.Basis1(), A.Basis1(), B.base().Basis1());
   // Iterate over the components in M, first index
   for (LinearAlgebra::const_iterator<OperatorComponent>::type I = iterate(M); I; ++I)
   {
      // second index in M
      for (LinearAlgebra::const_inner_iterator<OperatorComponent>::type J = iterate(I); J; ++J)
      {
         // Iterate over the irreducible components of M(I,J)
         for (SimpleRedOperator::const_iterator k = J->begin(); k != J->end(); ++k)
         {
            // *k is an irreducible operator.  Iterate over the components of this operator
            for (LinearAlgebra::const_iterator<SimpleOperator>::type R = iterate(*k); R; ++R)
            {
               for (LinearAlgebra::const_inner_iterator<SimpleOperator>::type 
                       S = iterate(R); S; ++S)
               {
                  add_triple_prod(Result[J.index1()], *S, 
                                  A[S.index1()], 
                                  E[J.index2()], 
                                  herm(B.base()[S.index2()]),
                                  k->TransformsAs(),
                                  M.Basis1()[J.index1()]);
               }
            }
         }
      }
   }
   return Result;
}

StateComponent
operator_prod(HermitianProxy<OperatorComponent> const& M,
              HermitianProxy<StateComponent> const& A, 
              StateComponent const& E,
              StateComponent const& B)
{
   DEBUG_PRECONDITION_EQUAL(M.base().LocalBasis2(), A.base().LocalBasis());
   DEBUG_PRECONDITION_EQUAL(M.base().LocalBasis1(), B.LocalBasis());
   DEBUG_PRECONDITION_EQUAL(M.base().Basis1(), E.LocalBasis());
   DEBUG_PRECONDITION_EQUAL(A.base().Basis1(), E.Basis1());
   DEBUG_PRECONDITION_EQUAL(E.Basis2(), B.Basis1());

   StateComponent Result(M.base().Basis2(), A.base().Basis2(), B.Basis2());
   // Iterate over the components in M, first index
   for (LinearAlgebra::const_iterator<OperatorComponent>::type I = iterate(M.base()); I; ++I)
   {
      // second index in M
      for (LinearAlgebra::const_inner_iterator<OperatorComponent>::type J = iterate(I); J; ++J)
      {
         // Iterate over the irreducible components of M(I,J)
         for (SimpleRedOperator::const_iterator k = J->begin(); k != J->end(); ++k)
         {
            // *k is an irreducible operator.  Iterate over the components of this operator
            for (LinearAlgebra::const_iterator<SimpleOperator>::type R = iterate(*k); R; ++R)
            {
               for (LinearAlgebra::const_inner_iterator<SimpleOperator>::type 
                       S = iterate(R); S; ++S)
               {
                  add_triple_prod(Result[J.index2()], herm(*S), 
                                  herm(A.base()[S.index1()]), 
                                  E[J.index1()], 
                                  B[S.index2()],
                                  k->TransformsAs(),
                                  M.base().Basis2()[J.index2()]);
               }
            }
         }
      }
   }
   return Result;
}

StateComponent
operator_prod_regular(OperatorComponent const& M,
		      StateComponent const& A, 
		      StateComponent const& F,
		      LinearAlgebra::HermitianProxy<StateComponent> const& B)
{
   return operator_prod(M, A, F, B);
}

StateComponent
operator_prod_regular(LinearAlgebra::HermitianProxy<OperatorComponent> const& M,
		      LinearAlgebra::HermitianProxy<StateComponent> const& A, 
		      StateComponent const& E,
		      StateComponent const& B)
{
   return operator_prod(M, A, E, B);
}

StateComponent
operator_prod_inner(OperatorComponent const& M,
                    StateComponent const& A, 
                    StateComponent const& E,
                    HermitianProxy<StateComponent> const& B)
{
   DEBUG_PRECONDITION_EQUAL(M.Basis1(), A.LocalBasis());
   DEBUG_PRECONDITION_EQUAL(M.Basis2(), B.base().LocalBasis());
   DEBUG_PRECONDITION_EQUAL(M.LocalBasis2(), E.LocalBasis());
   DEBUG_PRECONDITION_EQUAL(A.Basis2(), E.Basis1());
   DEBUG_PRECONDITION_EQUAL(E.Basis2(), B.base().Basis2());

   StateComponent Result(M.LocalBasis1(), A.Basis1(), B.base().Basis1());
   // Iterate over the components in M, first index
   for (LinearAlgebra::const_iterator<OperatorComponent>::type I = iterate(M); I; ++I)
   {
      // second index in M
      for (LinearAlgebra::const_inner_iterator<OperatorComponent>::type J = iterate(I); J; ++J)
      {
         // Iterate over the irreducible components of M(I,J)
         for (SimpleRedOperator::const_iterator k = J->begin(); k != J->end(); ++k)
         {
            // *k is an irreducible operator.  Iterate over the components of this operator
            for (LinearAlgebra::const_iterator<SimpleOperator>::type R = iterate(*k); R; ++R)
            {
               for (LinearAlgebra::const_inner_iterator<SimpleOperator>::type 
                       S = iterate(R); S; ++S)
               {
                  Result[S.index1()] += (*S) * triple_prod(A[J.index1()], 
                                                           E[S.index2()], 
                                                           herm(B.base()[J.index2()]),
                                                           k->TransformsAs(),
                                                           M.LocalBasis1()[S.index1()]);
               }
            }
         }
      }
   }
   return Result;
}

StateComponent
operator_prod_inner(HermitianProxy<OperatorComponent> const& M,
                    HermitianProxy<StateComponent> const& A, 
                    StateComponent const& E,
                    StateComponent const& B)
{
   DEBUG_PRECONDITION_EQUAL(M.base().Basis2(), A.base().LocalBasis());
   DEBUG_PRECONDITION_EQUAL(M.base().Basis1(), B.LocalBasis());
   DEBUG_PRECONDITION_EQUAL(M.base().LocalBasis1(), E.LocalBasis());
   DEBUG_PRECONDITION_EQUAL(A.base().Basis1(), E.Basis1());
   DEBUG_PRECONDITION_EQUAL(E.Basis2(), B.Basis1());

   StateComponent Result(M.base().LocalBasis2(), A.base().Basis2(), B.Basis2());
   // Iterate over the components in M, first index
   for (LinearAlgebra::const_iterator<OperatorComponent>::type I = iterate(M.base()); I; ++I)
   {
      // second index in M
      for (LinearAlgebra::const_inner_iterator<OperatorComponent>::type J = iterate(I); J; ++J)
      {
         // Iterate over the irreducible components of M(I,J)
         for (SimpleRedOperator::const_iterator k = J->begin(); k != J->end(); ++k)
         {
            // *k is an irreducible operator.  Iterate over the components of this operator
            for (LinearAlgebra::const_iterator<SimpleOperator>::type R = iterate(*k); R; ++R)
            {
               for (LinearAlgebra::const_inner_iterator<SimpleOperator>::type 
                       S = iterate(R); S; ++S)
               {
                  Result[S.index2()] += herm(*S) * triple_prod(herm(A.base()[J.index1()]), 
                                                               E[S.index1()], 
                                                               B[J.index2()],
                                                               k->TransformsAs(),
                                                               M.base().LocalBasis2()[S.index2()]);
               }
            }
         }
      }
   }
   return Result;
}

OperatorComponent
project_rows(OperatorComponent const& x, std::set<int> const& Rows)
{
   BasisList ProjectedBasis(x.Basis1().GetSymmetryList());
   std::vector<int> Map(x.Basis1().size(), -1);
   for (std::set<int>::const_iterator I = Rows.begin(); I != Rows.end(); ++I)
   {
      Map[*I] = ProjectedBasis.size();
      ProjectedBasis.push_back(x.Basis1()[*I]);
   }

   OperatorComponent Result(x.LocalBasis1(), x.LocalBasis2(), ProjectedBasis, x.Basis2());
   for (OperatorComponent::const_iterator I = iterate(x.data()); I; ++I)
   {
      for (OperatorComponent::const_inner_iterator J = iterate(I); J; ++J)
      {
         if (Map[J.index1()] >= 0)
         {
            Result(Map[J.index1()], J.index2()) = *J;
         }
      }
   }

   return Result;
}

// project onto the given columns of the component
OperatorComponent
project_columns(OperatorComponent const& x, std::set<int> const& Cols)
{
   BasisList ProjectedBasis(x.Basis2().GetSymmetryList());
   std::vector<int> Map(x.Basis2().size(), -1);
   for (std::set<int>::const_iterator I = Cols.begin(); I != Cols.end(); ++I)
   {
      Map[*I] = ProjectedBasis.size();
      ProjectedBasis.push_back(x.Basis2()[*I]);
   }

   OperatorComponent Result(x.LocalBasis1(), x.LocalBasis2(), x.Basis1(), ProjectedBasis);
   for (OperatorComponent::const_iterator I = iterate(x.data()); I; ++I)
   {
      for (OperatorComponent::const_inner_iterator J = iterate(I); J; ++J)
      {
         if (Map[J.index2()] >= 0)
         {
            Result(J.index1(), Map[J.index2()]) = *J;
         }
      }
   }

   return Result;
}

std::pair<OperatorComponent, OperatorComponent>
decompose_tensor_prod(SimpleOperator const& Op, 
                      ProductBasis<BasisList, BasisList> const& B1,
                      ProductBasis<BasisList, BasisList> const& B2)
{
   // the norm of the operator sets the scale for removing small components
   double NormScale = norm_frob(Op) / std::sqrt(double(Op.Basis1().size() * Op.Basis2().size()));

   // Firstly, decompose the tensor product of Op
   typedef std::map<Tensor::PartialProdIndex, std::complex<double> > PartialProdType;
   PartialProdType PartialProd = Tensor::decompose_tensor_prod(Op, B1, B2);

   // rearrange the matrix elements to do a singular value decomposition
   // To do this, we write the decomposed operator as a matrix indexed by
   // (qLeft,Left1,Left2), (qRight,Right1,Right2)
   // and do a SVD of it.  
   // Because we assume a scalar operator, qLeft is the adjoint of qRight,
   // and the resulting matrix is block diagonal with respect to qRight/qLeft.
   
   std::set<QuantumNumber> qUsed;
   for (PartialProdType::const_iterator I = PartialProd.begin(); I != PartialProd.end(); ++I)
   {
      qUsed.insert(I->first.qRight);
   }

   BasisList BondBasis(Op.GetSymmetryList()); // bond of the MPO
   std::vector<SimpleOperator> MatLeft, MatRight; // decompose into sum_i MatLeft[i] * MatRight[i]

   // for each quantum number (block diagonal of the partial transpose),
   // assemble the matrix indexed by (left1,left2) (right1,right2) and construct the SVD
   for (std::set<QuantumNumber>::const_iterator qUsedIter = qUsed.begin(); qUsedIter != qUsed.end(); ++qUsedIter)
   {
      // linearize the (left1,left2) and (right1,right2) indices
      LinearAlgebra::Matrix<int> LeftMap(B1.Left().size(), B2.Left().size(), -1);
      LinearAlgebra::Matrix<int> RightMap(B1.Right().size(), B2.Right().size(), -1);
      std::vector<std::pair<int, int> > LeftInverseMap, RightInverseMap;
      int LeftCount = 0, RightCount = 0;
      for (PartialProdType::const_iterator I = PartialProd.begin(); I != PartialProd.end(); ++I)
      {
         // only look at the components that transform properly
         if (I->first.qRight != *qUsedIter)
            continue;

         if (LeftMap(I->first.Left1, I->first.Left2) == -1)
         {
            LeftInverseMap.push_back(std::make_pair(I->first.Left1, I->first.Left2));
            LeftMap(I->first.Left1, I->first.Left2) = LeftCount++;
         }

         if (RightMap(I->first.Right1, I->first.Right2) == -1)
         {
            RightInverseMap.push_back(std::make_pair(I->first.Right1, I->first.Right2));
            RightMap(I->first.Right1, I->first.Right2) = RightCount++;
         }
      }

      // assemble the matrix
      LinearAlgebra::Matrix<std::complex<double> > Mat(LeftCount, RightCount, 0.0);
      for (PartialProdType::const_iterator I = PartialProd.begin(); I != PartialProd.end(); ++I)
      {
         // only look at the components that transform properly
         if (I->first.qRight != *qUsedIter)
            continue;

         Mat(LeftMap(I->first.Left1, I->first.Left2), RightMap(I->first.Right1, I->first.Right2)) = I->second;
      }

      // Perform the singular decomposition
      LinearAlgebra::Matrix<std::complex<double> > U, Vh;
      LinearAlgebra::Vector<double> D;
      SingularValueDecomposition(Mat, U, D, Vh);

      //CHECK(norm_frob(Mat - U*diagonal_matrix(D)*Vh) < 1E-10)(Mat)(U*diagonal_matrix(D)*Vh);

      // split into pairs of operators, keeping only the non-zero singular values
      for (unsigned k = 0; k < size(D); ++k)
      {
         if (D[k] > std::numeric_limits<double>::epsilon() * NormScale * 10)
         {
            double Coeff = std::sqrt(D[k]);  // distribute sqrt(D) to each operator
            SimpleOperator L(B1.Left(), B2.Left(), adjoint(*qUsedIter));
            for (unsigned x = 0; x < size1(U); ++x)
            {
               if (norm_frob(U(x,k)) > std::numeric_limits<double>::epsilon() * 10)
               {
                  L(LeftInverseMap[x].first, LeftInverseMap[x].second) = U(x,k) * Coeff;
               }
            }

            SimpleOperator R(B1.Right(), B2.Right(), *qUsedIter);
            for (unsigned x = 0; x < size2(Vh); ++x)
            {
               if (norm_frob(Vh(k,x)) > std::numeric_limits<double>::epsilon() * 10)
               {
                  R(RightInverseMap[x].first, RightInverseMap[x].second) = Vh(k,x) * Coeff;
               }
            }

            BondBasis.push_back(*qUsedIter);  // add the entry to the bond basis for this operator
            MatLeft.push_back(L);
            MatRight.push_back(R);
         }
      }
   }

   // Now we have the final operator decomposition and bond basis, construct the MPO's.

   OperatorComponent MA(B1.Left(), B2.Left(), make_vacuum_basis(Op.GetSymmetryList()), BondBasis);
   OperatorComponent MB(B1.Right(), B2.Right(), BondBasis, make_vacuum_basis(Op.GetSymmetryList()));
   for (unsigned i = 0; i < BondBasis.size(); ++i)
   {
      MA(0,i) = MatLeft[i];
      MB(i,0) = MatRight[i];
   }

   // test
#if !defined(NDEBUG)
   SimpleOperator OpTest = local_tensor_prod(MA, MB)(0,0).scalar();
   CHECK(norm_frob(OpTest - Op) < 1E-13 * NormScale)(Op)(OpTest)(norm_frob(OpTest-Op));
#endif

   return std::pair<OperatorComponent, OperatorComponent>(MA, MB);
}

std::pair<OperatorComponent, OperatorComponent>
decompose_tensor_prod(SimpleOperator const& Op, 
                      ProductBasis<BasisList, BasisList> const& B)
{
   return decompose_tensor_prod(Op, B, B);
}

std::pair<OperatorComponent, OperatorComponent>
decompose_tensor_prod(SimpleOperator const& Op, 
                      BasisList const& BasisA, BasisList const& BasisB)
{
   ProductBasis<BasisList, BasisList> B(BasisA, BasisB);
   return decompose_tensor_prod(Op, B, B);
}

OperatorComponent
RotateToOperatorRightBoundary(StateComponent const& x)
{
   IrredTensor<LinearAlgebra::Matrix<double>, VectorBasis, BasisList> Splitter1 = SplitBasis(x.Basis1());
   IrredTensor<LinearAlgebra::Matrix<double>, VectorBasis, BasisList> Splitter2 = SplitBasis(x.Basis2());

   OperatorComponent Result(Splitter1.Basis2(), Splitter2.Basis2(), x.LocalBasis(), make_vacuum_basis(x.GetSymmetryList()));

   for (unsigned i = 0; i < x.size(); ++i)
   {
      Result(i, 0) = map_1x1_operator(triple_prod(herm(Splitter1), x[i], Splitter2));
   }

   return Result;
}

OperatorComponent
RotateToOperatorLeftBoundary(StateComponent const& x)
{
   IrredTensor<LinearAlgebra::Matrix<double>, VectorBasis, BasisList> Splitter1 = SplitBasis(x.Basis1());
   IrredTensor<LinearAlgebra::Matrix<double>, VectorBasis, BasisList> Splitter2 = SplitBasis(x.Basis2());

   OperatorComponent Result(adjoint(Splitter1.Basis2()), adjoint(Splitter2.Basis2()), 
                            make_vacuum_basis(x.GetSymmetryList()), x.LocalBasis());

   for (unsigned i = 0; i < x.size(); ++i)
   {
      Result(0, i) = flip_conj(map_1x1_operator(triple_prod(herm(Splitter1), x[i], Splitter2)));
   }

   return Result;
}

OperatorComponent
conj(OperatorComponent const& x)
{
   OperatorComponent Result(x);
   for (OperatorComponent::iterator I = iterate(Result.data()); I; ++I)
   {
      for (OperatorComponent::inner_iterator J = iterate(I); J; ++J)
      {
         *J = conj(*J);
      }
   }
   return Result;
}

OperatorComponent
flip_conj(OperatorComponent const& x)
{
   OperatorComponent Result(adjoint(x.LocalBasis1()), adjoint(x.LocalBasis2()), adjoint(x.Basis1()), adjoint(x.Basis2()));
   for (OperatorComponent::const_iterator I = iterate(x.data()); I; ++I)
   {
      for (OperatorComponent::const_inner_iterator J = iterate(I); J; ++J)
      {
         TRACE(*J)(flip_conj(*J));
         Result(J.index1(), J.index2()) = flip_conj(*J);
      }
   }
   return Result;
}

void
update_mask_basis1(std::vector<int>& Mask1, SimpleOperator const& Op, std::vector<int> const& Mask2)
{
   for (const_iterator<SimpleOperator>::type I = iterate(Op); I; ++I)
   {
      for (const_inner_iterator<SimpleOperator>::type J = iterate(I); J; ++J)
      {
         if (Mask2[J.index2()])
            Mask1[J.index1()] = true;
      }
   }
}

void
update_mask_basis1(std::vector<int>& Mask1, SimpleRedOperator const& Op, std::vector<int> const& Mask2)
{
   for (SimpleRedOperator::const_iterator I = Op.begin(); I != Op.end(); ++I)
   {
      update_mask_basis1(Mask1, *I, Mask2);
   }
}

void
update_mask_basis1(std::vector<int>& Mask1, OperatorComponent const& Op, std::vector<int> const& Mask2)
{
   for (OperatorComponent::const_iterator I = iterate(Op.data()); I; ++I)
   {
      for (OperatorComponent::const_inner_iterator J = iterate(I); J; ++J)
      {
         update_mask_basis1(Mask1, *J, Mask2);
      }
   }
}

void
update_mask_basis2(std::vector<int> const& Mask1, SimpleOperator const& Op, std::vector<int>& Mask2)
{
   for (const_iterator<SimpleOperator>::type I = iterate(Op); I; ++I)
   {
      for (const_inner_iterator<SimpleOperator>::type J = iterate(I); J; ++J)
      {
         if (Mask1[J.index1()])
            Mask2[J.index2()] = true;
      }
   }
}

void
update_mask_basis2(std::vector<int> const& Mask1, SimpleRedOperator const& Op, std::vector<int>& Mask2)
{
   for (SimpleRedOperator::const_iterator I = Op.begin(); I != Op.end(); ++I)
   {
      update_mask_basis2(Mask1, *I, Mask2);
   }
}

void
update_mask_basis2(std::vector<int> const& Mask1, OperatorComponent const& Op, std::vector<int>& Mask2)
{
   for (OperatorComponent::const_iterator I = iterate(Op.data()); I; ++I)
   {
      for (OperatorComponent::const_inner_iterator J = iterate(I); J; ++J)
      {
         update_mask_basis2(Mask1, *J, Mask2);
      }
   }
}

MatrixOperator ExpandBasis1Used(StateComponent& A, OperatorComponent const& Op)
{
   DEBUG_CHECK_EQUAL(A.LocalBasis(), Op.LocalBasis1());
   DEBUG_CHECK_EQUAL(A.LocalBasis(), Op.LocalBasis2());
   std::vector<int> Mask(A.size(), 0);
   for (unsigned i = 0; i < A.size(); ++i)
   {
      Mask[i] = (norm_frob(A[i]) > std::numeric_limits<double>::epsilon());
   }

   std::vector<int> Mask1 = Mask;
   update_mask_basis1(Mask1, Op, Mask);

   return ExpandBasis1Used(A, Mask1);
}

MatrixOperator ExpandBasis2Used(StateComponent& A, OperatorComponent const& Op)
{
   DEBUG_CHECK_EQUAL(A.LocalBasis(), Op.LocalBasis1());
   DEBUG_CHECK_EQUAL(A.LocalBasis(), Op.LocalBasis2());
   std::vector<int> Mask(A.size(), 0);
   for (unsigned i = 0; i < A.size(); ++i)
   {
      Mask[i] = (norm_frob(A[i]) > std::numeric_limits<double>::epsilon());
   }

   std::vector<int> Mask2 = Mask;
   update_mask_basis2(Mask, Op, Mask2);

   return ExpandBasis2Used(A, Mask2);
}
