// -*- C++ -*- $Id$

#include "operator_component.h"
#include "tensor/tensorproduct.h"
#include "linearalgebra/eigen.h"
#include "tensor/regularize.h"
#include "linearalgebra/matrix_utility.h"

//using namespace LinearAlgebra;
using LinearAlgebra::operator*;
using LinearAlgebra::adjoint;

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

double
norm_frob_sq(OperatorComponent const& x)
{
   return norm_frob_sq(local_inner_prod(herm(x), x));
}

void
OperatorComponent::check_structure() const
{
   for (const_iterator I = iterate(*this); I; ++I)
   {
      for (const_inner_iterator J = iterate(I); J; ++J)
      {
         for (SimpleRedOperator::const_iterator q = J->begin(); q != J->end(); ++q)
         {
            CHECK(is_transform_target(this->Basis2()[J.index2()], q->TransformsAs(), this->Basis1()[J.index1()]))
               (this->Basis1())(J.index1())(this->Basis2())(J.index2())(q->TransformsAs());
         }
      }
   }
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
operator-(OperatorComponent const& A, OperatorComponent const& B)
{
   DEBUG_CHECK_EQUAL(A.LocalBasis1(), B.LocalBasis1());
   DEBUG_CHECK_EQUAL(A.LocalBasis2(), B.LocalBasis2());
   DEBUG_CHECK_EQUAL(A.Basis1(), B.Basis1());
   DEBUG_CHECK_EQUAL(A.Basis2(), B.Basis2());
   OperatorComponent Result(A.LocalBasis1(), A.LocalBasis2(), A.Basis1(), A.Basis2());
   Result.data() = A.data() - B.data();
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
   Result.debug_check_structure();
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
   Result.debug_check_structure();
   return Result;
}
 
OperatorComponent 
tensor_row_sum(OperatorComponent const& A, 
	       OperatorComponent const& B)
{
   if (A.is_null())
      return B;
   if (B.is_null())
      return A;
   return tensor_row_sum(A, B, SumBasis<BasisList>(A.Basis2(), B.Basis2()));
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
   Result.debug_check_structure();
   return Result;
}

OperatorComponent 
tensor_col_sum(OperatorComponent const& A, 
	       OperatorComponent const& B)
{
   if (A.is_null())
      return B;
   if (B.is_null())
      return A;
   return tensor_col_sum(A, B, SumBasis<BasisList>(A.Basis1(), B.Basis1()));
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

SimpleOperator
local_inner_prod(HermitianProxy<OperatorComponent> const& A, OperatorComponent const& B)
{
   DEBUG_PRECONDITION_EQUAL(A.base().Basis1(), B.Basis1());

   SimpleOperator Result(A.base().Basis2(), B.Basis2());

   OperatorComponent::const_iterator AI = iterate(A.base());
   OperatorComponent::const_iterator BI = iterate(B);
   while (AI)
   {
      for (OperatorComponent::const_inner_iterator AJ = iterate(AI); AJ; ++AJ)
      {
         for (OperatorComponent::const_inner_iterator BJ = iterate(BI); BJ; ++BJ)
         {
            Result(AJ.index2(), BJ.index2()) = inner_prod(*AJ, *BJ);
         }
      }
      ++AI;
      ++BI;
   }
   Result.debug_check_structure();
   return Result;
}


SimpleOperator
local_inner_prod(OperatorComponent const& A, HermitianProxy<OperatorComponent> const& B)
{
   DEBUG_PRECONDITION_EQUAL(A.Basis2(), B.base().Basis2());

   SimpleOperator Result(A.Basis1(), B.base().Basis1());

   for (OperatorComponent::const_iterator AI = iterate(A); AI; ++AI)
   {
      for (OperatorComponent::const_iterator BI = iterate(B.base()); BI; ++BI)
      {
         OperatorComponent::const_inner_iterator AJ = iterate(AI);
         OperatorComponent::const_inner_iterator BJ = iterate(BI);

         while (AJ && BJ)
         {
            if (AJ.index2() == BJ.index2())
            {
               Result(AJ.index1(), BJ.index1()) += conj(inner_prod(*BJ, *AJ));
	       ++AJ;
	       ++BJ;
            }
            else if (AJ.index2() < BJ.index2())
            {
               ++AJ;
            }
            else
            {
               ++BJ;
            }
         }
      }
   }
   Result.debug_check_structure();
   return Result;
}

OperatorComponent
local_adjoint(OperatorComponent const& A)
{
   OperatorComponent Result(A.LocalBasis2(), A.LocalBasis1(), adjoint(A.Basis1()), adjoint(A.Basis2()));
   for (OperatorComponent::const_iterator I = iterate(A); I; ++I)
   {
      for (OperatorComponent::const_inner_iterator J = iterate(I); J; ++J)
      {
         Result(J.index1(), J.index2()) = adjoint(*J);
      }
   }
   Result.debug_check_structure();
   return Result;
}

OperatorComponent aux_tensor_prod(OperatorComponent const& ML, OperatorComponent const& MR)
{
   DEBUG_PRECONDITION_EQUAL(ML.LocalBasis2(), MR.LocalBasis1());
   ProductBasis<BasisList, BasisList> PBasis1(ML.Basis1(), MR.Basis1());
   ProductBasis<BasisList, BasisList> PBasis2(ML.Basis2(), MR.Basis2());

   OperatorComponent Result(ML.LocalBasis1(), MR.LocalBasis2(), PBasis1.Basis(), PBasis2.Basis());

   for (OperatorComponent::const_iterator iL = iterate(ML); iL; ++iL)
   {
      for (OperatorComponent::const_iterator iR = iterate(MR); iR; ++iR)
      {

         for (OperatorComponent::const_inner_iterator jL = iterate(iL); jL; ++jL)
         {
            for (OperatorComponent::const_inner_iterator jR = iterate(iR); jR; ++jR)
            {

               ProductBasis<BasisList, BasisList>::const_iterator TiIter, 
                  TiEnd = PBasis1.end(jL.index1(), jR.index1());
               for (TiIter = PBasis1.begin(jL.index1(), jR.index1()); TiIter != TiEnd; ++TiIter)
               {
		  ProductBasis<BasisList, BasisList>::const_iterator TjIter, 
		     TjEnd = PBasis2.end(jL.index2(), jR.index2());
		  for (TjIter = PBasis2.begin(jL.index2(), jR.index2()); TjIter != TjEnd; 
		       ++TjIter)
                  {
                     for (SimpleRedOperator::const_iterator RedL = jL->begin(); RedL != jL->end(); ++RedL)
                     {
                        for (SimpleRedOperator::const_iterator RedR = jR->begin(); RedR != jR->end(); ++RedR)
                        {
                           QuantumNumbers::QuantumNumberList FinalQ = transform_targets(RedL->TransformsAs(), RedR->TransformsAs());
                           for (QuantumNumbers::QuantumNumberList::const_iterator q = FinalQ.begin(); q != FinalQ.end(); ++q)
                           {
                              if (is_transform_target(PBasis2[*TjIter], *q, PBasis1[*TiIter]))
                              {
                                 double Coeff = tensor_coefficient(PBasis1, PBasis2, 
                                                                   RedL->TransformsAs(), RedR->TransformsAs(), *q,
                                                                   jL.index1(), jR.index1(), *TiIter, 
                                                                   jL.index2(), jR.index2(), *TjIter);
                                 if (LinearAlgebra::norm_2(Coeff) > 1E-14)
                                 {
                                    Result(*TiIter, *TjIter) += Coeff * prod(*RedL, *RedR, *q);
                                 }
                              }
                           }
                        }
                     }
		  }
	       }
	    }
	 }
      }
   }
   Result.debug_check_structure();
   return Result;
}

#if 0
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
   Result.debug_check_structure();
   return Result;
}

SimpleOperator TruncateBasis1(OperatorComponent& A)
{
   // We want to work from the last row to the first, so that we preserve the last row exactly.
   // We don't have to worry about accidentally eliminating the first row of a triangular MPO,
   // since if it is triangular then the last column cannot be parallel to anything else

   // iterate through the rows of A and remove any that are parallel with a previous.
   // We cannot safely remove rows that are merely linearly dependent on the previous
   // rows without potential trouble with catastrophic cancellation.

   // if the operator is trivially zero, then return early
   if (A.Basis1().size() == 0)
   {
      return SimpleOperator(A.Basis1(), A.Basis1(), QuantumNumber(A.GetSymmetryList()));
   }
   else if (A.Basis2().size() == 0) 
   {
      SimpleOperator Result(A.Basis1(), A.Basis2(), QuantumNumber(A.GetSymmetryList()));
      A = OperatorComponent(A.LocalBasis1(), A.LocalBasis2(), A.Basis2(), A.Basis2());
      return Result;
   }

   int Size = A.Basis1().size();

   // Calculate the matrix of overlaps of each row
   SimpleOperator Overlaps = local_inner_prod(A, herm(A));

   // This is the transform that truncates the rows of A.
   // row[i] = NewRows[i].second * A(NewRows[i].second, all)
   // The flag value -1 indicates that the row is not needed
   std::vector<std::pair<int, std::complex<double> > > NewRows(Size, std::make_pair(-1, 1.0));
   std::set<int> KeepStates; // the states that we are going to keep
   for (int i = Size-1; i >= 0; --i)
   {
      NewRows[i] = std::make_pair(i, 1.0);
      double imat = Overlaps(i,i).real();
      // if the row is zero, we can eliminate it completely.
      // Because we've normalized everything, then it is either 1 or ~epsilon here.
      if (imat == 0)
         continue;    // skip this row

      bool Parallel = false;  // is vector i parallel to some other vector?
      // loop to find out if row i is parallel to row j
      for (int j = A.Basis1().size()-1; j > i; --j)
      {
         // can only have non-zero determinant if the quantum numbers match
         if (A.Basis1()[i] != A.Basis1()[j])
            continue;

         // If row j was removed from the basis, then we are not parallel to it
         if (NewRows[j].first == -1)
            continue;

         double const OverlapEpsilon = 1E-14;

         double jmat = Overlaps(j,j).real();
         std::complex<double> ijmat = Overlaps(i,j);
         // are rows i and j parallel?
         if ((imat*jmat - LinearAlgebra::norm_frob_sq(ijmat)) / (imat*imat + jmat*jmat) < OverlapEpsilon)
         {
            NewRows[i] = std::make_pair(j, ijmat / jmat);
            // corner case: i is parallel to j, but we might have previously determined that
            // j is parallel to some other row k.  Does this ever happen in practice?
            if (NewRows[j].first != j)
            {
               WARNING("parallel row vectors have a non-transitive equivalence")(i)(j)(NewRows[j].first);
               DEBUG_TRACE(A);
               while (NewRows[i].first != i && NewRows[NewRows[i].first].first != NewRows[i].first)
               {
                  NewRows[i] = std::make_pair(NewRows[NewRows[i].first].first,
                                              NewRows[NewRows[i].first].second * NewRows[i].second);
               }
            }
            Parallel = true;
            break;  // we have a parallel vector, break out of the loop
         }
      }
      if (!Parallel)
         KeepStates.insert(i);
   }

   //TRACE_IF(KeepStates.size() != A.Basis1().size());

   // construct the reduced basis
   BasisList NewBasis(A.GetSymmetryList());
   std::map<unsigned, unsigned> BasisMapping;  // maps the kept states into the new basis
   for (std::set<int>::const_iterator I = KeepStates.begin(); I != KeepStates.end(); ++I)
   {
     BasisMapping[*I] = NewBasis.size();
     NewBasis.push_back(A.Basis1()[*I]);
   }

   // construct the projector from the old onto the new basis
   SimpleOperator Trunc(NewBasis, A.Basis1(), QuantumNumber(A.GetSymmetryList()));
   for (std::set<int>::const_iterator I = KeepStates.begin(); I != KeepStates.end(); ++I)
   {
      Trunc(BasisMapping[*I], *I) = 1.0;
   }

   // construct the regularizer to give the full basis from the reduced basis
   SimpleOperator Reg(A.Basis1(), NewBasis, QuantumNumber(A.GetSymmetryList()));
   for (unsigned i = 0; i < NewRows.size(); ++i)
   {
      if (NewRows[i].first != -1)
         Reg(i, BasisMapping[NewRows[i].first]) = NewRows[i].second;
   }

   OperatorComponent tA = prod(Trunc, A); // the result

   //#if !defined(NDEBUG)
   // verify that prod(Reg, tA) is the same as A.  
   OperatorComponent ACheck = prod(Reg, tA);
   CHECK(norm_frob(A - ACheck) < 1E-10)(tA)(ACheck)(Trunc)(Reg)(Overlaps);
   //#endif

   A = tA;
   return Reg;
}

SimpleOperator TruncateBasis2(OperatorComponent& A)
{
   // We want to work from the first column to last, so that we preserve the first column exactly.
   // We don't have to worry about accidentally eliminating the last column of a triangular MPO,
   // since if it is triangular then the last column cannot be parallel to anything else

   // if the operator is trivially zero, then return early
   if (A.Basis2().size() == 0)
   {
      return SimpleOperator(A.Basis2(), A.Basis2(), QuantumNumber(A.GetSymmetryList()));
   }
   else if (A.Basis1().size() == 0)
   {
      SimpleOperator Result(A.Basis1(), A.Basis2(), QuantumNumber(A.GetSymmetryList()));
      A = OperatorComponent(A.LocalBasis1(), A.LocalBasis2(), A.Basis1(), A.Basis1());
      return Result;
   }

   SimpleOperator Overlaps = local_inner_prod(herm(A), A);

   // This is the transform that truncates the columns of A.
   // row[i] = NewCols[i].second * A(all, NewCols[i].second)
   std::vector<std::pair<int, std::complex<double> > > NewCols;
   std::set<int> KeepStates; // the states that we are going to keep
   for (int i = 0; i < int(A.Basis2().size()); ++i)
   {
      NewCols.push_back(std::make_pair(i, 1.0));
      double imat = Overlaps(i,i).real();
      // if the row is zero, we can eliminate it completely
      if (imat == 0)
      {
         NewCols.back().first = -1;  // special value, indicates the row is not needed
         continue;
      }
      bool Parallel = false;  // is vector i parallel to some other vector?
      // loop to find out if row i is parallel to row j
      for (int j = 0; j < i; ++j)
      {
         // can only have non-zero determinant if the quantum numbers match
         if (A.Basis2()[i] != A.Basis2()[j])
            continue;

         // If column j was removed from the basis, then we are not parallel to it
         if (NewCols[j].first == -1)
            continue;

         double const OverlapEpsilon = 1E-14;

         double jmat = Overlaps(j,j).real();
         std::complex<double> ijmat = Overlaps(j,i);
         // are rows i and j parallel?
         if ((imat*jmat - LinearAlgebra::norm_frob_sq(ijmat)) / (imat*imat + jmat*jmat) < OverlapEpsilon)
         {
            NewCols[i] = std::make_pair(j, ijmat / jmat);
            // corner case: i is parallel to j, but we might have previously determined that
            // j is parallel to some other column k.  Does this ever happen in practice?
            if (NewCols[j].first != j)
            {
               WARNING("parallel column vectors have a non-transitive equivalence")(i)(j)(NewCols[j].first);
               DEBUG_TRACE(A);
               while (NewCols[i].first != i && NewCols[NewCols[i].first].first != NewCols[i].first)
               {
                  NewCols[i] = std::make_pair(NewCols[NewCols[i].first].first,
                                              NewCols[NewCols[i].first].second * NewCols[i].second);
               }
            }
            Parallel = true;
            break;  // we have a parallel vector, break out of the loop
         }
      }
      if (!Parallel)
         KeepStates.insert(i);
   }

   //TRACE_IF(KeepStates.size() != A.Basis1().size());

   // construct the reduced basis
   BasisList NewBasis(A.GetSymmetryList());
   std::map<unsigned, unsigned> BasisMapping;  // maps the kept states into the new basis
   for (std::set<int>::const_iterator I = KeepStates.begin(); I != KeepStates.end(); ++I)
   {
     BasisMapping[*I] = NewBasis.size();
     NewBasis.push_back(A.Basis2()[*I]);
   }

   // construct the projector from the old onto the new basis
   SimpleOperator Trunc(A.Basis2(), NewBasis, QuantumNumber(A.GetSymmetryList()));
   for (std::set<int>::const_iterator I = KeepStates.begin(); I != KeepStates.end(); ++I)
   {
      Trunc(*I, BasisMapping[*I]) = 1.0;
   }

   // construct the regularizer to give the full basis from the reduced basis
   SimpleOperator Reg(NewBasis, A.Basis2(), QuantumNumber(A.GetSymmetryList()));
   for (unsigned i = 0; i < NewCols.size(); ++i)
   {
      if (NewCols[i].first != -1)
         Reg(BasisMapping[NewCols[i].first], i) = NewCols[i].second;
   }

   // adjust for the normalizer
   //   Reg = Reg * InverseNormalizer;

   OperatorComponent tA = prod(A, Trunc); // the result

#if !defined(NDEBUG)
   // verify that prod(tA, Reg) is the same as A.  
   OperatorComponent ACheck = prod(tA, Reg);
   CHECK(norm_frob(A - ACheck) < 1E-10)(tA)(ACheck)(Trunc)(Reg)(Overlaps);
#endif

   A = tA;
   return Reg;
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
   DEBUG_PRECONDITION_EQUAL(M.base().LocalBasis2(), A.base().LocalBasis());
   DEBUG_PRECONDITION_EQUAL(M.base().LocalBasis1(), B.LocalBasis());
   DEBUG_PRECONDITION_EQUAL(M.base().Basis1(), E.LocalBasis());
   DEBUG_PRECONDITION_EQUAL(A.base().Basis1(), E.Basis1());
   DEBUG_PRECONDITION_EQUAL(E.Basis2(), B.Basis1());

   DEBUG_PRECONDITION_EQUAL(A.base().Basis2(), B.Basis2());
   DEBUG_PRECONDITION_EQUAL(E.Basis1(), E.Basis2());

   StateComponent Result(M.base().Basis2(), A.base().Basis2(), B.Basis2());

   // First component is identity
   Result.front() = MatrixOperator::make_identity(B.Basis2());
   for (unsigned alpha = 1; alpha < M.base().Basis2().size(); ++alpha)
   {
      Result[alpha] = MatrixOperator(B.Basis2(), B.Basis2(), M.base().Basis2()[alpha]);
      // only need to iterate to alpha, since upper triangular

      // alphap = 0 part, in this case we assume E[0] = identity so simplifies further
      if (M.base().iterate_at(0, alpha))
      {
         // TODO: we could further optimize this if M(0,alpha) is proportional to identity,
         // but this only happens if the MPO is non-optimsl?
         Result[alpha] += operator_prod(herm(M.base()(0,alpha)), A, B);
      }
      for (unsigned alphap = 1; alphap < alpha; ++alphap)
      {
         if (M.base().iterate_at(alphap, alpha))
         {
            Result[alpha] += operator_prod(herm(M.base()(alphap, alpha)), A, E[alphap], B, M.base().Basis2()[alpha]);
         }
      }
   }
   return Result;
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

      // This sets the scale of the singular values.  Any singular values smaller than
      // KeepThreshold are removed.
      double KeepThreshold = std::numeric_limits<double>::epsilon() * LinearAlgebra::max(D) * 10;

      // split into pairs of operators, keeping only the non-zero singular values
      for (unsigned k = 0; k < size(D); ++k)
      {
         if (D[k] > KeepThreshold)
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
   CHECK(norm_frob(OpTest - Op) < 1E-13 * norm_frob(OpTest))(Op)(OpTest)(norm_frob(OpTest-Op));
#endif

   return std::pair<OperatorComponent, OperatorComponent>(MA, MB);
}

std::pair<OperatorComponent, OperatorComponent>
decompose_tensor_prod(OperatorComponent const& Op, 
                      ProductBasis<BasisList, BasisList> const& B1,
                      ProductBasis<BasisList, BasisList> const& B2)
{
   BasisList Vacuum(Op.GetSymmetryList());
   OperatorComponent Left(B1.Left(), B2.Left(), Op.Basis1(), Vacuum);
   OperatorComponent Right(B1.Right(), B2.Right(), Vacuum, Op.Basis2());
   // we treat this as a sum and separately decompose each part
   for (OperatorComponent::const_iterator I = iterate(Op); I; ++I)
   {
      for (OperatorComponent::const_inner_iterator J = iterate(I); J; ++J)
      {
	 // Make 1xN and Nx1 OperatorComponent's for *J
	 OperatorComponent L, R;
	 for (SimpleRedOperator::const_iterator RI = J->begin(); RI != J->end(); ++RI)
	 {
	    OperatorComponent l, r;
	    boost::tie(l, r) = decompose_tensor_prod(*RI, B1, B2);
	    L = tensor_row_sum(L, l);
	    R = tensor_col_sum(R, r);
	 }
	 // move (L,R) into the (i,j) component
	 OperatorComponent TL(B1.Left(), B2.Left(), Op.Basis1(), L.Basis2());
	 for (OperatorComponent::const_iterator T1 = iterate(L); T1; ++T1)
	 {
	    for (OperatorComponent::const_inner_iterator T2 = iterate(T1); T2; ++T2)
	    {
	       TL(J.index1(), T2.index2()) = *T2;
	    }
	 }
	 OperatorComponent TR(B1.Right(), B2.Right(), R.Basis1(), Op.Basis2());
	 for (OperatorComponent::const_iterator T1 = iterate(R); T1; ++T1)
	 {
	    for (OperatorComponent::const_inner_iterator T2 = iterate(T1); T2; ++T2)
	    {
	       TR(T2.index1(), J.index2()) = *T2;
	    }
	 }
	 // Now we can sum to the result
	 Left = tensor_row_sum(Left, TL);
	 Right = tensor_col_sum(Right, TR);
      }
   }
   // This is likely to leave the middle bond basis quite big, so optimize it.
   // (We maybe need to iterate this function?)
   SimpleOperator T = TruncateBasis1(Right);
   Left = Left * T;
   T = TruncateBasis2(Left);
   Right = T * Right;
   return std::make_pair(Left, Right);
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
         DEBUG_TRACE(*J)(flip_conj(*J));
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

SimpleOperator CollapseBasis(BasisList const& b)
{
   std::set<QuantumNumber> QN(b.begin(), b.end());
   BasisList NewB(b.GetSymmetryList(), QN.begin(), QN.end());
   SimpleOperator C(NewB, b);
   for (unsigned j = 0; j < b.size(); ++j)
   {
      unsigned i = std::find(NewB.begin(), NewB.end(), b[j]) - NewB.begin();
      C(i,j) = 1.0;
   }
   return C;
}

SimpleOperator ProjectBasis(BasisList const& b, QuantumNumbers::QuantumNumber const& q)
{
   BasisList NewB(q);
   SimpleOperator Result(NewB, b);
   for (unsigned j = 0; j < b.size(); ++j)
   {
      if (b[j] == q)
         Result(0, j) = 1.0;
   }
   return Result;
}
