// -*- C++ -*- $Id$

#include "operator_component.h"
#include "tensor/tensorproduct.h"
#include "linearalgebra/eigen.h"
#include "tensor/regularize.h"
#include "linearalgebra/matrix_utility.h"
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

// Epsilon for truncation overlaps.
// The scale of Epsilon here is the square of the machine precision
// This needs some deeper analysis - the problem is we get numerically near parallel
// vectors although one of them is small and mostly noise, so we then end up with
// some spurious large matrix elements.
double const TruncateOverlapEpsilon = 1E-18;

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
               (this->Basis1())(J.index1())(this->Basis2())(J.index2())(q->TransformsAs())(*this);
         }
      }
   }
}

SimpleOperator
swap_gate(BasisList const& B1, BasisList const& B2, 
	  ProductBasis<BasisList, BasisList> const& Basis_21,
	  ProductBasis<BasisList, BasisList> const& Basis_12)
{
   QuantumNumbers::QuantumNumber Ident(B1.GetSymmetryList());

   SimpleOperator Result(Basis_21.Basis(), Basis_12.Basis(), Ident);
   for (unsigned i = 0; i < Basis_12.size(); ++i)
   {
      std::pair<int,int> x = Basis_12.rmap(i);

      // Find the corresponding element in Basis1 with swapped indices
      for (ProductBasis<BasisList>::const_iterator I = Basis_21.begin(x.second, x.first); I != Basis_21.end(x.second, x.first); ++I)
      {
	 if (Basis_21[*I] == Basis_12[i])
	 {
	    // This phase factor works for the xxx spin chain
	    Result(*I, i) = conj_phase(Basis_12[i], adjoint(Basis_12.Left()[x.first]), Basis_12.Right()[x.second]);
	 }
      }
   }
   //TRACE(Result)(Result.Basis1())(Result.Basis2());
   return Result;
}

SimpleOperator
swap_gate_fermion(BasisList const& B1, LinearAlgebra::Vector<double> const& Parity1, 
		  BasisList const& B2, LinearAlgebra::Vector<double> const& Parity2,
		  ProductBasis<BasisList, BasisList> const& Basis_21,
		  ProductBasis<BasisList, BasisList> const& Basis_12)
{
   DEBUG_CHECK_EQUAL(Parity1.size(), B1.size());
   DEBUG_CHECK_EQUAL(Parity2.size(), B2.size());
   QuantumNumbers::QuantumNumber Ident(B1.GetSymmetryList());

   SimpleOperator Result(Basis_21.Basis(), Basis_12.Basis(), Ident);
   for (unsigned i = 0; i < Basis_12.size(); ++i)
   {
      std::pair<int,int> x = Basis_12.rmap(i);

      // Find the corresponding element in Basis1 with swapped indices
      for (ProductBasis<BasisList>::const_iterator I = Basis_21.begin(x.second, x.first); I != Basis_21.end(x.second, x.first); ++I)
      {
	 if (Basis_21[*I] == Basis_12[i])
	 {
	    // This phase factor works for the xxx spin chain
	    Result(*I, i) = conj_phase(Basis_12[i], adjoint(Basis_12.Left()[x.first]), Basis_12.Right()[x.second])
	       * Parity1[x.first] * Parity2[x.second];
	 }
      }
   }
   //TRACE(Result)(Result.Basis1())(Result.Basis2());
   return Result;
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

void print_structure(OperatorComponent const& Op, std::ostream& out, double UnityEpsilon)
{
   for (unsigned i = 0; i < Op.size1(); ++i)
   {
      out << '[';
      for (unsigned j = 0; j < Op.size2(); ++j)
      {
	 if (Op.iterate_at(i,j))
	 {
	    SimpleRedOperator X = Op(i,j);
	    if (X.size() > 1)
	       out << 'x';       // some compound operator
	    else if (!is_scalar(X.begin()->TransformsAs()))
	    {
	       out << 'v';       // a non-scalar
	    }
	    else
	    {
	       SimpleOperator Y = X.scalar();
	    
	       std::complex<double> x = PropIdent(Y);
	       if (x == 0.0)
	       {
		  std::complex<double> x = PropIdent(scalar_prod(herm(Y),Y), UnityEpsilon);
		  if (norm_frob(x-1.0) < 1E-12)
		     out << 'U';      // a unitary
		  else
		     out << 's';      // a generic scalar
	       }
	       else if (norm_frob(x-1.0) < 1E-12)
	       {
		  out << 'I';         // the identity
	       }
	       else
		  out << 'i';         // something proportional to the identity
	    }
	 }
	 else
	    out << ' ';
      }
      out << "]\n";
   }
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

   for (const_iterator<MatType>::type AI = iterate(A); AI; ++AI)
   {
      for (const_inner_iterator<MatType>::type AJ = iterate(AI); AJ; ++AJ)
      {
	 const_iterator<MatType>::type BI = iterate(B);
	 BI += AJ.index2();
	 for (const_inner_iterator<MatType>::type BJ = iterate(BI); BJ; ++BJ)
	 {
	    for (SimpleRedOperator::const_iterator AComp = AJ->begin(); AComp != AJ->end(); ++AComp)
	    {
	       for (SimpleRedOperator::const_iterator BComp = BJ->begin(); BComp != BJ->end(); ++BComp)
	       {
		  QuantumNumbers::QuantumNumberList qL = transform_targets(BComp->TransformsAs(), AComp->TransformsAs());
		  for (QuantumNumbers::QuantumNumberList::const_iterator qI = qL.begin(); qI != qL.end(); ++qI)
		  {
		     // Ensure that the final operator is a possible target
		     if (is_transform_target(B.Basis2()[BJ.index2()], *qI, A.Basis1()[AJ.index1()]))
		     {
			double Coeff = product_coefficient(AComp->TransformsAs(), BComp->TransformsAs(), *qI,
							   A.Basis1()[AJ.index1()], B.Basis2()[BJ.index2()], B.Basis1()[BJ.index1()]);
			if (LinearAlgebra::norm_2(Coeff) > 1E-14)
			   Result(AJ.index1(), BJ.index2()) += Coeff * tensor_prod(*AComp, *BComp, *qI);
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

SimpleOperator
local_inner_prod(HermitianProxy<OperatorComponent> const& A, OperatorComponent const& B)
{
   DEBUG_PRECONDITION_EQUAL(A.base().Basis1(), B.Basis1());

   // The coupling coefficient for the auxiliary basis is a simple rescaling, 
   // same as for the tensor ScalarProd operation

   SimpleOperator Result(A.base().Basis2(), B.Basis2());

   // AI/BI point at the same quantum number, and is the inner summation of the product
   OperatorComponent::const_iterator AI = iterate(A.base());
   OperatorComponent::const_iterator BI = iterate(B);
   while (AI)
   {
      double InnerDegree = degree(B.Basis1()[BI.index()]);

      for (OperatorComponent::const_inner_iterator AJ = iterate(AI); AJ; ++AJ)
      {
         for (OperatorComponent::const_inner_iterator BJ = iterate(BI); BJ; ++BJ)
         {
	    if (A.base().Basis2()[AJ.index2()] == B.Basis2()[BJ.index2()])
	       Result(AJ.index2(), BJ.index2()) += InnerDegree / degree(Result.Basis1()[AJ.index2()]) * inner_prod(*AJ, *BJ);
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

   // FIXME: is there another coupling coefficient here for the product?
   // That we are taking a scalar product simplifies it a bit, but there should still be a coefficient.

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
	       if (A.Basis1()[AJ.index1()] == B.base().Basis1()[BJ.index1()])
		  Result(AJ.index1(), BJ.index1()) += adjoint(inner_prod(*BJ, *AJ));
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

SimpleOperator
local_inner_tensor_prod(HermitianProxy<OperatorComponent> const& A, OperatorComponent const& B)
{
   DEBUG_PRECONDITION_EQUAL(A.base().Basis1(), B.Basis1());

   ProductBasis<BasisList, BasisList> PBasis1(adjoint(A.base().Basis1()), B.Basis1());
   ProductBasis<BasisList, BasisList> PBasis2(adjoint(A.base().Basis2()), B.Basis2());

   QuantumNumbers::QuantumNumber Ident(B.GetSymmetryList());

   SimpleOperator Result(PBasis1.Basis(), PBasis2.Basis(), Ident);

   for (OperatorComponent::const_iterator iL = iterate(A.base()); iL; ++iL)
   {
      for (OperatorComponent::const_iterator iR = iterate(B); iR; ++iR)
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
			   // we are taking the inner product, so we require the quantum numbers match
			   if (RedL->TransformsAs() != RedR->TransformsAs())
			      continue;
			   
			   if (!is_transform_target(PBasis2[*TjIter], Ident, PBasis1[*TiIter]))
			      continue;

			   double Coeff = tensor_coefficient(PBasis1, PBasis2, 
							     adjoint(RedL->TransformsAs()), RedR->TransformsAs(), Ident,
							     jL.index1(), jR.index1(), *TiIter, 
							     jL.index2(), jR.index2(), *TjIter);
			   if (LinearAlgebra::norm_2(Coeff) > 1E-14)
			   {
			      Result(*TiIter, *TjIter) += Coeff * inner_prod(*RedL, *RedR);
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

   // A norm for the overlaps matrix
   double Scale = trace(Overlaps).real() / Overlaps.Basis1().size();

   // This is the transform that truncates the rows of A.
   // row[i] = NewRows[i].second * A(NewRows[i].second, all)
   // The flag value -1 indicates that the row is not needed
   std::vector<std::pair<int, std::complex<double> > > NewRows(Size, std::make_pair(-1, 1.0));
   std::set<int> KeepStates; // the states that we are going to keep
   for (int i = Size-1; i >= 0; --i)
   {
      double imat = Overlaps(i,i).real();
      // if the row is numerically zero, we can eliminate it completely.
      if (imat <= Scale * TruncateOverlapEpsilon*TruncateOverlapEpsilon)
         continue;    // skip this row

      NewRows[i] = std::make_pair(i, 1.0);

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

         double jmat = Overlaps(j,j).real();
         std::complex<double> ijmat = Overlaps(i,j);
         // are rows i and j parallel?
         if ((imat*jmat - LinearAlgebra::norm_frob_sq(ijmat)) / (imat*imat + jmat*jmat) < TruncateOverlapEpsilon)
         {
            NewRows[i] = std::make_pair(j, conj(ijmat) / jmat);
            // corner case: i is parallel to j, but we might have previously determined that
            // j is parallel to some other row k.  Does this ever happen in practice?
            if (NewRows[j].first != j)
            {
               WARNING("parallel row vectors have a non-transitive equivalence")(i)(j)(NewRows[j].first)
		  (Overlaps(j,i))(Overlaps(NewRows[j].first,i))(Overlaps(NewRows[j].first,j))(imat)(jmat);
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
   CHECK(norm_frob_sq(A - ACheck) <= Scale*TruncateOverlapEpsilon)(norm_frob_sq(A - ACheck))(A)(ACheck)(A-ACheck)(Trunc)(Reg)(Overlaps);
   //#endif

   A = tA;
   return Reg;
}

SimpleOperator TruncateBasis2(OperatorComponent& A)
{
   // We want to work from the first column to last, so that we preserve the first column exactly.
   // For a TriangularMPO the first column will contain the identity.
   // We don't have to worry about accidentally eliminating the last column of a triangular MPO,
   // since if it is triangular then the last column cannot be parallel to anything else.

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

   // A norm for the overlaps matrix
   double Scale = trace(Overlaps).real() / Overlaps.Basis1().size();

   // This is the transform that truncates the columns of A.
   // row[i] = NewCols[i].second * A(all, NewCols[i].second)
   std::vector<std::pair<int, std::complex<double> > > NewCols;
   std::set<int> KeepStates; // the states that we are going to keep
   for (int i = 0; i < int(A.Basis2().size()); ++i)
   {
      NewCols.push_back(std::make_pair(i, 1.0));
      double imat = Overlaps(i,i).real();
      // if the row is zero, we can eliminate it completely
      if (imat <= Scale * TruncateOverlapEpsilon*TruncateOverlapEpsilon)
      {
         NewCols.back().first = -1;  // special value, indicates the row is not needed
         continue;
      }
      bool Parallel = false;  // is vector i parallel to some other vector?
      // loop to find out if row i is parallel to some previous column j
      for (int j = 0; j < i; ++j)
      {
         // can only have non-zero determinant if the quantum numbers match
         if (A.Basis2()[i] != A.Basis2()[j])
            continue;

         // If column j was removed from the basis, then we are not parallel to it
         if (NewCols[j].first == -1)
            continue;

         double jmat = Overlaps(j,j).real();
         std::complex<double> ijmat = Overlaps(j,i);
         // are rows i and j parallel?
         if ((imat*jmat - LinearAlgebra::norm_frob_sq(ijmat)) / (imat*imat + jmat*jmat) < TruncateOverlapEpsilon)
         {
            NewCols[i] = std::make_pair(j, ijmat / jmat);
            // corner case: i is parallel to j, but we might have previously determined that
            // j is parallel to some other column k, whereas we found that column i is NOT
	    // parallel to k.
            if (NewCols[j].first != j)
            {
	       int k = (NewCols[j].first);
               WARNING("parallel column vectors have a non-transitive equivalence")(i)(j)(NewCols[j].first)
			(Overlaps(j,i))(Overlaps(k,i))(Overlaps(k,j))(Overlaps(k,k))(imat)(jmat);
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

   //#if !defined(NDEBUG)
   // verify that prod(tA, Reg) is the same as A.  
   OperatorComponent ACheck = prod(tA, Reg);
   CHECK(norm_frob_sq(A - ACheck) <= (Scale*TruncateOverlapEpsilon))(norm_frob_sq(A - ACheck))(A)(ACheck)(A-ACheck)(Trunc)(Reg)(Overlaps);
   //#endif

   A = tA;
   return Reg;
}

SimpleOperator TruncateBasis1MkII(OperatorComponent& A, double Epsilon)
{
   // We want to work from the last row to the first, so that we preserve the last row exacrly.
   // For a TriangularMPO the last row will contain the identity.
   // We don't have to worry about accidentally eliminating the first row of a triangular MPO,
   // since if it is triangular then the last column cannot be parallel to anything else.

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

   // A norm for the overlaps matrix
   double Scale = norm_frob_sq(A) / (A.Basis1().total_degree() * A.Basis2().total_degree());

   // make a dense matrix
   LinearAlgebra::Matrix<SimpleRedOperator> M = A.data();

   // these are stored in the reverse order
   std::vector<LinearAlgebra::Vector<SimpleRedOperator> > Rows;
   std::vector<LinearAlgebra::Vector<std::complex<double> > > T;

   std::vector<double> RowNormSq;

   std::vector<QuantumNumbers::QuantumNumber> NewBasis1Q;

   double Normalization = A.LocalBasis2().total_degree();

   int r = M.size1()-1;
   while (r >= 0)
   {
      LinearAlgebra::Vector<SimpleRedOperator> NextRow = M(r, LinearAlgebra::all);
      
      // orthogonalize against the following rows
      for (unsigned i = 0; i < Rows.size(); ++i)
      {
	 // check that the quantum numbers agree - for non-abelian quantum numbers it might
	 // have non-zero overlap but be in a different symmetry sector
	 if (NewBasis1Q[i] == A.Basis1()[r])
	 {
	    std::complex<double> x = inner_prod(Rows[i], NextRow) / RowNormSq[i];
	    T[i][r] += x;
	    if (norm_frob(x) > 0) // Epsilon*Epsilon)
	       NextRow = NextRow - LinearAlgebra::Vector<SimpleRedOperator>(x * Rows[i]);
	 }
      }
      
      // DGKS step
      for (unsigned i = 0; i < Rows.size(); ++i)
      {
	 // check that the quantum numbers agree - for non-abelian quantum numbers it might
	 // have non-zero overlap but be in a different symmetry sector
	 if (NewBasis1Q[i] == A.Basis1()[r])
	 {
	    std::complex<double> x = inner_prod(Rows[i], NextRow) / RowNormSq[i];
	    T[i][r] += x;
	    if (norm_frob(x) > 0) // Epsilon*Epsilon)
	       NextRow = NextRow - LinearAlgebra::Vector<SimpleRedOperator>(x * Rows[i]);
	 }
      }

      // check the norm of column c
      double NextRowNormSq = norm_frob_sq(NextRow);
      if (NextRowNormSq > Epsilon*Epsilon)
      {
	 // we have a column
	 Rows.push_back((std::sqrt(Normalization / NextRowNormSq)) * NextRow);
	 T.push_back(LinearAlgebra::Vector<std::complex<double> >(M.size1(), 0.0));
	 T[Rows.size()-1][r] = std::sqrt(NextRowNormSq / Normalization);
	 NewBasis1Q.push_back(A.Basis1()[r]);
	 RowNormSq.push_back(Normalization);  // since we've already normalized it
      }

      --r;
   }

   BasisList NewBasis1(NewBasis1Q.rbegin(), NewBasis1Q.rend());

   // Convert back to OperatorComponent and SimpleOperator
   OperatorComponent ANew(A.LocalBasis1(), A.LocalBasis2(), NewBasis1, A.Basis2());
   for (unsigned r = 0; r < NewBasis1.size(); ++r)
   {
      for (unsigned c = 0; c < A.Basis2().size(); ++c)
      {
	 if (norm_frob(Rows[r][c]) > 0)
	 {
	    ANew.data()(NewBasis1.size()-r-1,c) = Rows[r][c];
	    ANew.data()(NewBasis1.size()-r-1,c).trim();
	 }
      }
   }

   SimpleOperator Trunc(A.Basis1(), NewBasis1);
   for (unsigned r = 0; r < A.Basis1().size(); ++r)
   {
      for (unsigned d = 0; d < NewBasis1.size(); ++d)
      {
	 if (norm_frob(T[NewBasis1.size()-d-1][r]) > 0)
	 {
	    Trunc.data()(r,d) = T[NewBasis1.size()-d-1][r];
	 }
      }
   }

   OperatorComponent ACheck = prod(Trunc, ANew);
   CHECK(norm_frob(A - ACheck) <= Scale*TruncateOverlapEpsilon);
   
   ANew.check_structure();
   Trunc.check_structure();

   A = ANew;
   return Trunc;
}

SimpleOperator TruncateBasis2MkII(OperatorComponent& A, double Epsilon)
{
   // We want to work from the first column to last, so that we preserve the first column exactly.
   // For a TriangularMPO the first column will contain the identity.
   // We don't have to worry about accidentally eliminating the last column of a triangular MPO,
   // since if it is triangular then the last column cannot be parallel to anything else.

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

   // A norm for the overlaps matrix
   double Scale = norm_frob_sq(A) / (A.Basis1().total_degree() * A.Basis2().total_degree());

   // make a dense matrix
   LinearAlgebra::Matrix<SimpleRedOperator> M = A.data();

   std::vector<LinearAlgebra::Vector<SimpleRedOperator> > Columns;
   std::vector<LinearAlgebra::Vector<std::complex<double> > > T;

   std::vector<double> ColNormSq;

   BasisList NewBasis2(A.Basis2().GetSymmetryList());

   double Normalization = A.LocalBasis1().total_degree();

   int c = 0;
   while (c < int(M.size2()))
   {
      //      LinearAlgebra::Vector<std::complex<double> > v(M.size2(), 0.0);
      
      LinearAlgebra::Vector<SimpleRedOperator> NextCol = M(LinearAlgebra::all, c);
      
      // orthogonalize against the previous rows
      for (unsigned i = 0; i < Columns.size(); ++i)
      {
	 // check that the quantum numbers agree - for non-abelian quantum numbers it might
	 // have non-zero overlap but be in a different symmetry sector
	 if (NewBasis2[i] == A.Basis2()[c])
	 {
	    std::complex<double> x = inner_prod(Columns[i], NextCol) /  ColNormSq[i];
	    T[i][c] += x;
	    if (norm_frob(x) > 0)
	       //Epsilon*Epsilon)
	       NextCol = NextCol - LinearAlgebra::Vector<SimpleRedOperator>(x * Columns[i]);
	 }
      }
      
      // DGKS step
      for (unsigned i = 0; i < Columns.size(); ++i)
      {
	 // check that the quantum numbers agree - for non-abelian quantum numbers it might
	 // have non-zero overlap but be in a different symmetry sector
	 if (NewBasis2[i] == A.Basis2()[c])
	 {
	    std::complex<double> x = inner_prod(Columns[i], NextCol) /  ColNormSq[i];
	    T[i][c] += x;
	    if (norm_frob(x) > 0)
	       //Epsilon*Epsilon)
	       NextCol = NextCol - LinearAlgebra::Vector<SimpleRedOperator>(x * Columns[i]);
	 }
      }

      // check the norm of column c
      double NextColNormSq = norm_frob_sq(NextCol);
      if (NextColNormSq > Epsilon*Epsilon)
      {
	 // we have a column
	 Columns.push_back((std::sqrt(Normalization / NextColNormSq)) * NextCol);
	 T.push_back(LinearAlgebra::Vector<std::complex<double> >(M.size2(), 0.0));
	 T[Columns.size()-1][c] = std::sqrt(NextColNormSq / Normalization);
	 NewBasis2.push_back(A.Basis2()[c]);
	 ColNormSq.push_back(Normalization);  // since we've already normalized it
      }

      ++c;
   }

   // Convert back to OperatorComponent and SimpleOperator
   OperatorComponent ANew(A.LocalBasis1(), A.LocalBasis2(), A.Basis1(), NewBasis2);
   for (unsigned r = 0; r < A.Basis1().size(); ++r)
   {
      for (unsigned c = 0; c < NewBasis2.size(); ++c)
      {
	 if (norm_frob(Columns[c][r]) > 0)
	 {
	    ANew.data()(r,c) = Columns[c][r];
	    ANew.data()(r,c).trim();
	 }
      }
   }

   SimpleOperator Trunc(NewBasis2, A.Basis2());
   for (unsigned c = 0; c < NewBasis2.size(); ++c)
   {
      for (unsigned d = 0; d < A.Basis2().size(); ++d)
      {
	 if (norm_frob(T[c][d]) > 0)
	 {
	    Trunc.data()(c,d) = T[c][d];
	 }
      }
   }

   OperatorComponent ACheck = prod(ANew, Trunc);
   CHECK(norm_frob(A - ACheck) <= Scale*TruncateOverlapEpsilon);

   ANew.check_structure();
   Trunc.check_structure();

   A = ANew;
   return Trunc;
}

StateComponent
operator_prod(OperatorComponent const& M,
              StateComponent const& A, 
              StateComponent const& E,
              HermitianProxy<StateComponent> const& B)
{
   PRECONDITION_EQUAL(M.LocalBasis1(), A.LocalBasis());
   PRECONDITION_EQUAL(M.LocalBasis2(), B.base().LocalBasis());
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
   PRECONDITION_EQUAL(M.base().LocalBasis2(), A.base().LocalBasis());
   PRECONDITION_EQUAL(M.base().LocalBasis1(), B.LocalBasis());
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
   PRECONDITION_EQUAL(M.base().LocalBasis2(), A.base().LocalBasis());
   PRECONDITION_EQUAL(M.base().LocalBasis1(), B.LocalBasis());
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

std::vector<std::pair<SimpleRedOperator, SimpleRedOperator> >
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
   
   BasisList BondBasis(Op.GetSymmetryList()); // bond of the MPO
   std::vector<SimpleRedOperator> MatLeft, MatRight; // decompose into sum_i MatLeft[i] * MatRight[i]

   // linearize the tuples (qLeft,Left1,Left2) and (qRight,Right1,Right2) into integers.
   // We also keep the inverse map, so we can quickly convert the integer back into a tuple
   typedef boost::tuple<QuantumNumbers::QuantumNumber, int, int> PartialProdIndex;
   std::map<PartialProdIndex, int> LeftIndex, RightIndex;
   std::vector<PartialProdIndex> LeftReverse, RightReverse;
   for (PartialProdType::const_iterator I = PartialProd.begin(); I != PartialProd.end(); ++I)
   { 
      PartialProdIndex Left(I->first.qLeft, I->first.Left1, I->first.Left2);
      // If Left doesn't exist yet, then add it as the next successive integer (given by LeftIndex.size())
      if (LeftIndex.find(Left) == LeftIndex.end())
      {
	 LeftIndex[Left] = LeftReverse.size();
	 LeftReverse.push_back(Left);
      }
	 
      // same for the right index
      PartialProdIndex Right(I->first.qRight, I->first.Right1, I->first.Right2);
      // If Right doesn't exist yet, then add it as the next successive integer (given by RightIndex.size())
      if (RightIndex.find(Right) == RightIndex.end())
      {
	 RightIndex[Right] = RightReverse.size();
	 RightReverse.push_back(Right);
      }
   }

   // Now assemble the matrix
   LinearAlgebra::Matrix<std::complex<double> > Mat(LeftIndex.size(), RightIndex.size(), 0.0);
   for (PartialProdType::const_iterator I = PartialProd.begin(); I != PartialProd.end(); ++I)
   {
      PartialProdIndex Left(I->first.qLeft, I->first.Left1, I->first.Left2);
      PartialProdIndex Right(I->first.qRight, I->first.Right1, I->first.Right2);
      Mat(LeftIndex[Left], RightIndex[Right]) = I->second;
   }
   
   // Perform the singular decomposition
   LinearAlgebra::Matrix<std::complex<double> > U, Vh;
   LinearAlgebra::Vector<double> D;
   SingularValueDecomposition(Mat, U, D, Vh);

   // This sets the scale of the singular values.  Any singular values smaller than
   // KeepThreshold are removed.
   double KeepThreshold = std::numeric_limits<double>::epsilon() * LinearAlgebra::max(D) * 10;
      
   // split into pairs of operators, keeping only the non-zero singular values
   std::vector<std::pair<SimpleRedOperator, SimpleRedOperator> > Result;
   for (unsigned k = 0; k < size(D); ++k)
   {
      DEBUG_TRACE(D[k]);
      if (D[k] > KeepThreshold)
      {
	 double Coeff = std::sqrt(D[k]);  // distribute sqrt(D) to each operator
	 SimpleRedOperator L(B1.Left(), B2.Left());
	 for (unsigned x = 0; x < size1(U); ++x)
	 {
	    if (norm_frob(U(x,k)) > std::numeric_limits<double>::epsilon() * 10)
	    {
	       L.project(LeftReverse[x].get<0>())
		  (LeftReverse[x].get<1>(), LeftReverse[x].get<2>()) = U(x,k) * Coeff;
	    }
	 }
	 
	 SimpleRedOperator R(B1.Right(), B2.Right());
	 for (unsigned x = 0; x < size2(Vh); ++x)
	 {
	    if (norm_frob(Vh(k,x)) > std::numeric_limits<double>::epsilon() * 10)
	    {
	       R.project(RightReverse[x].get<0>())
		  (RightReverse[x].get<1>(), RightReverse[x].get<2>()) = Vh(k,x) * Coeff;
	    }
	 }

	 Result.push_back(std::make_pair(L, R));
      }
   }

#if !defined(NDEBUG)
   SimpleOperator OpTest;
   for (unsigned i = 0; i < Result.size(); ++i)
   {
      // We really need an IrredProd for reducible tensors, it would be more efficient here
      OpTest += project(Result[i].first * Result[i].second, Op.TransformsAs());
   }
   CHECK(norm_frob(OpTest - Op) < 1E-13 * norm_frob(OpTest))(Op)(OpTest)(norm_frob(OpTest-Op))
      (MatLeft)(MatRight);
#endif

   return Result;
}

std::pair<OperatorComponent, OperatorComponent>
decompose_local_tensor_prod(OperatorComponent const& Op, 
			    ProductBasis<BasisList, BasisList> const& B1,
			    ProductBasis<BasisList, BasisList> const& B2)
{
   // A scale factor for determining whether to keep numerically small values
   double ScaleFactor = norm_frob_sq(Op) / std::sqrt(Op.Basis1().size() * Op.Basis2().size() * Op.LocalBasis1().size() * Op.LocalBasis2().size());

   typedef std::map<Tensor::PartialProdIndex, std::complex<double> > PartialProdType;

   // ComponentIndex is an index to a particular matrix element of the MPO transforming as q,
   // between local indices local1,local2 and auxiliary index aux
   // ComponentIndex(aux, local1, local2, q)
   typedef boost::tuple<int, int, int, QuantumNumbers::QuantumNumber> ComponentIndex;

   // The full index stores the ComponentIndex for the left and right operators
   typedef std::pair<ComponentIndex, ComponentIndex> FullIndex;

   // Represent a decomposed tensor product for a single internal quantum number
   typedef std::map<FullIndex, std::complex<double> > DecomposedMPOType;

   // and finally sum over possible internal quantum numbers
   typedef std::map<QuantumNumbers::QuantumNumber, DecomposedMPOType> DecomposedByQuantumType;
   DecomposedByQuantumType DecomposedOperator;

   // Mapping from the ComponentIndex to a linear integer index, for the left and right operators
   std::map<QuantumNumbers::QuantumNumber, std::map<ComponentIndex, int> > LeftMapping, RightMapping;
   // Reverse mapping from integer index to ComponentIndex for left and right operators
   std::map<QuantumNumbers::QuantumNumber, std::vector<ComponentIndex> > LeftRevMapping, RightRevMapping;

   for (OperatorComponent::const_iterator I = iterate(Op); I; ++I)
   {
      int LeftAux1 = I.index();
      for (OperatorComponent::const_inner_iterator J = iterate(I); J; ++J)
      {
	 int RightAux2 = J.index2();
	 for (SimpleRedOperator::const_iterator RI = J->begin(); RI != J->end(); ++RI)
	 {
	    // Decompose the tensor product component
	    PartialProdType PartialProd = Tensor::decompose_tensor_prod(*RI, B1, B2);

	    // decompose the product in the auxiliary basis
	    // This is the operator <q' | AB(k) | q > where
	    // q' = Basis1[J.index1()]
	    // q = Basis2[J.index2()]
	    // k = RI->TransformsAs()
	    // decomposed into < q' | A(k1) | q'' > < q'' | B(k2) | q > = sum_k <q' | AB(k) | q >
	    // with 
	    // k1 = PartialProd.qLeft
	    // k2 = PartialProd.qRight
	    // and we sum over q'', being all possible internal quantum numbers,
	    // and this determines internal basis.

	    // For each possible q'', we assemble the matrix, into DecomposedOperator[q'']
	    for (PartialProdType::const_iterator P = PartialProd.begin(); P != PartialProd.end(); ++P)
	    {
	       QuantumNumbers::QuantumNumberList QL 
		  = transform_targets(Op.Basis2()[RightAux2], P->first.qRight);

	       for (QuantumNumbers::QuantumNumberList::const_iterator Qpp = QL.begin(); Qpp != QL.end(); ++Qpp)
	       {
		  // Skip quantum numbers that are not possible (structural zero in the coupling coefficient)
		  if (!is_transform_target(*Qpp,  P->first.qLeft, Op.Basis1()[LeftAux1]))
		     continue;

		  // Qpp is the quantum number of the inner basis

		  double Coeff = inverse_product_coefficient(P->first.qLeft, P->first.qRight, RI->TransformsAs(),
							     Op.Basis1()[J.index1()], Op.Basis2()[J.index2()], *Qpp);
		  if (LinearAlgebra::norm_2(Coeff) > 1E-14)
		  {
		     ComponentIndex LeftIndex 
			= boost::make_tuple(LeftAux1, P->first.Left1, P->first.Left2,
					    P->first.qLeft);
		     if (LeftMapping[*Qpp].find(LeftIndex) == LeftMapping[*Qpp].end())
		     {
			LeftMapping[*Qpp][LeftIndex] = LeftRevMapping[*Qpp].size();
			LeftRevMapping[*Qpp].push_back(LeftIndex);
		     }

		     ComponentIndex RightIndex
			= boost::make_tuple(RightAux2, P->first.Right1, P->first.Right2,
					    P->first.qRight);
		     if (RightMapping[*Qpp].find(RightIndex) == RightMapping[*Qpp].end())
		     {
			RightMapping[*Qpp][RightIndex] = RightRevMapping[*Qpp].size();
			RightRevMapping[*Qpp].push_back(RightIndex);
		     }

		     DecomposedOperator[*Qpp][std::make_pair(LeftIndex, RightIndex)]
			= Coeff * P->second;
		  }
	       }
	    }
	 }
      }
   }

   // Now we loop over the quantum numbers q'', assemble each matrix and do the SVD
   // We can't assemble this directly into MPOs since we need to know the internal basis first.
   // So instead we keep a list of kept singular vectors

   typedef std::map<ComponentIndex, std::complex<double> > SingularVectorType;
   typedef std::vector<SingularVectorType> SingularVectorListType;
   std::map<QuantumNumbers::QuantumNumber, SingularVectorListType> LeftSingularVectors, RightSingularVectors;

   BasisList Inner(Op.GetSymmetryList()); // we can assemble the inner basis at the same time

   // For deciding which singular values to keep, we need a scale factor, so we need the largest singular value.
   // First pass, calculate the largest singular value
   double LargestSingular = 0.0;
   for (DecomposedByQuantumType::const_iterator Q = DecomposedOperator.begin(); 
	Q != DecomposedOperator.end(); ++Q)
   {
      QuantumNumbers::QuantumNumber Qpp = Q->first;

      // assemble the matrix
      LinearAlgebra::Matrix<std::complex<double> > Mat(LeftMapping[Qpp].size(), RightMapping[Qpp].size(), 0.0);
      for (DecomposedMPOType::const_iterator I = Q->second.begin(); I != Q->second.end(); ++I)
      {
	 Mat(LeftMapping[Qpp][I->first.first], RightMapping[Qpp][I->first.second]) = I->second;
      }

      // Do the SVD
      LinearAlgebra::Matrix<std::complex<double> > U, Vh;
      LinearAlgebra::Vector<double> D;
      SingularValueDecomposition(Mat, U, D, Vh);

      double m = LinearAlgebra::max(D);
      if (m > LargestSingular)
	 LargestSingular = m;
   }

   //TRACE(LargestSingular);

   // This sets the scale of the singular values.  Any singular values smaller than
   // KeepThreshold are removed.
   double KeepThreshold = std::numeric_limits<double>::epsilon() * LargestSingular * 100;

   // Second pass, collate the kept singular values
   for (DecomposedByQuantumType::const_iterator Q = DecomposedOperator.begin(); 
	Q != DecomposedOperator.end(); ++Q)
   {
      QuantumNumbers::QuantumNumber Qpp = Q->first;

      // assemble the matrix
      LinearAlgebra::Matrix<std::complex<double> > Mat(LeftMapping[Qpp].size(), RightMapping[Qpp].size(), 0.0);
      for (DecomposedMPOType::const_iterator I = Q->second.begin(); I != Q->second.end(); ++I)
      {
	 Mat(LeftMapping[Qpp][I->first.first], RightMapping[Qpp][I->first.second]) = I->second;
      }

      // Do the SVD
      LinearAlgebra::Matrix<std::complex<double> > U, Vh;
      LinearAlgebra::Vector<double> D;
      SingularValueDecomposition(Mat, U, D, Vh);

      // Assemble the singular vectors
      for (unsigned k = 0; k < size(D); ++k)
      {
	 // Skip over vectors that don't reach the threshold
	 //	 if (D[k] < KeepThreshold)
	 //	    continue;

	 double Coeff = std::sqrt(D[k]);  // distribute sqrt(D) to each operator

	 SingularVectorType LeftVec;
	 for (unsigned x = 0; x < size1(U); ++x)
	 {
	    // We might have some components that are in forbidden quantum number sectors - these should be small.
	    // In the auxiliary space, we should have a matrix element that transforms as q,
	    // with q1 = aux, q2 = Qpp
	    // In the local space, we have a matrix element that transforms as q,
	    // with q1 = local1, q2 = local2,
	    // where
	    // LeftRevMapping[Qpp][x] is a ComponentIndex(aux, local1, local2, q)

	    if (!is_transform_target(Qpp, LeftRevMapping[Qpp][x].get<3>(), Op.Basis1()[LeftRevMapping[Qpp][x].get<0>()])
		|| !is_transform_target(B2.Left()[LeftRevMapping[Qpp][x].get<2>()], 
					LeftRevMapping[Qpp][x].get<3>(), 
					B1.Left()[LeftRevMapping[Qpp][x].get<1>()]))
	    {
	       TRACE("Ignoring forbidden off-diagonal matrix element")
		  (U(x,k))(Qpp)(LeftRevMapping[Qpp][x].get<3>())(Op.Basis1()[LeftRevMapping[Qpp][x].get<0>()])
		  (B2.Left()[LeftRevMapping[Qpp][x].get<2>()])(Qpp)(B1.Left()[LeftRevMapping[Qpp][x].get<1>()]);
	    }
	    else if (norm_frob(U(x,k)) > std::numeric_limits<double>::epsilon() * 10)
	    {
	       LeftVec[LeftRevMapping[Qpp][x]] = Coeff*U(x,k);
	    }
	 }

	 SingularVectorType RightVec;
	 for (unsigned x = 0; x < size2(Vh); ++x)
	 {
	    if (norm_frob(Vh(k,x)) > std::numeric_limits<double>::epsilon() * 10)
	    {
	       RightVec[RightRevMapping[Qpp][x]] = Coeff*Vh(k,x);
	    }
	 }

	 LeftSingularVectors[Qpp].push_back(LeftVec);
	 RightSingularVectors[Qpp].push_back(RightVec);

	 // Add the quantum number to the inner basis
	 Inner.push_back(Qpp);
      }
   }

   // Now the inner basis is complete we can construct the OperatorComponents
   OperatorComponent MA(B1.Left(), B2.Left(), Op.Basis1(), Inner);
   OperatorComponent MB(B1.Right(), B2.Right(), Inner, Op.Basis2());

   int b = 0; // linear index into the Inner basis
   for (DecomposedByQuantumType::const_iterator Q = DecomposedOperator.begin(); 
	Q != DecomposedOperator.end(); ++Q)
   {
      QuantumNumbers::QuantumNumber Qpp = Q->first;
      for (unsigned i = 0; i < LeftSingularVectors[Qpp].size(); ++i)
      {
	 SingularVectorType LeftVec = LeftSingularVectors[Qpp][i];
	 for (SingularVectorType::const_iterator I = LeftVec.begin(); I != LeftVec.end(); ++I)
	 {
	    project(MA(I->first.get<0>(), b), I->first.get<3>())
	       (I->first.get<1>(), I->first.get<2>())
	       = I->second;
	 }

	 SingularVectorType RightVec = RightSingularVectors[Qpp][i];
	 for (SingularVectorType::const_iterator I = RightVec.begin(); I != RightVec.end(); ++I)
	 {
	    project(MB(b, I->first.get<0>()), I->first.get<3>())
	       (I->first.get<1>(), I->first.get<2>())
	       = I->second;
	 }

	 ++b;
      }
   }

   MA.debug_check_structure();
   MB.debug_check_structure();

#if 0
   // This is likely to leave the middle bond basis quite big, so optimize it.
   // (We maybe need to iterate this function?)
   SimpleOperator T = TruncateBasis1(Right);
   Left = Left * T;
   T = TruncateBasis2(Left);
   Right = T * Right;
#endif

#if !defined(NDEBUG)
   OperatorComponent OpTest = local_tensor_prod(MA, MB);
   CHECK(norm_frob(OpTest - Op) < 1E-13 * norm_frob(OpTest))(Op)(OpTest)(norm_frob(OpTest-Op))
      (MA)(MB);
#endif

   return std::make_pair(MA, MB);
}

std::pair<OperatorComponent, OperatorComponent>
decompose_local_tensor_prod(SimpleOperator const& Op, 
			    ProductBasis<BasisList, BasisList> const& B1,
			    ProductBasis<BasisList, BasisList> const& B2)
{
   OperatorComponent OpC(Op.Basis1(), Op.Basis2(),
			 make_single_basis(Op.TransformsAs()), make_vacuum_basis(Op.GetSymmetryList()));
   project(OpC(0,0), Op.TransformsAs()) = Op;
   return decompose_local_tensor_prod(OpC, B1, B2);
}

std::pair<OperatorComponent, OperatorComponent>
decompose_local_tensor_prod(SimpleOperator const& Op, 
			    ProductBasis<BasisList, BasisList> const& B)
{
   return decompose_local_tensor_prod(Op, B, B);
}

std::pair<OperatorComponent, OperatorComponent>
decompose_local_tensor_prod(SimpleOperator const& Op, 
			    BasisList const& BasisA, BasisList const& BasisB)
{
   ProductBasis<BasisList, BasisList> B(BasisA, BasisB);
   return decompose_local_tensor_prod(Op, B, B);
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

std::complex<double> PropIdent(SimpleOperator const& X, double UnityEpsilon)
{
   DEBUG_PRECONDITION_EQUAL(X.Basis1(), X.Basis2());
   SimpleOperator Ident = SimpleOperator::make_identity(X.Basis1());
   std::complex<double> x = inner_prod(Ident, X) / double(X.Basis1().total_degree());
   if (norm_frob_sq(X-x*Ident) > UnityEpsilon*UnityEpsilon)
      x = 0.0;
   return x;
}

