// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mpo/operator_component.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Reseach publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

#include "operator_component.h"
#include "tensor/tensorproduct.h"
#include "tensor/regularize.h"
#include <tuple>

#include "common/environment.h"

// Epsilon for truncation overlaps.
// The scale of Epsilon here is the square of the machine precision
// This needs some deeper analysis - the problem is we get numerically near parallel
// vectors although one of them is small and mostly noise, so we then end up with
// some spurious large matrix elements.
double const TruncateOverlapEpsilon = getenv_or_default("MP_TOE", 1E-15);

double const QREpsilon = getenv_or_default("MP_QR_EPS", 1E-13);

double const RemoveRowEpsilon = getenv_or_default("MP_RRE", 1E-40);

#if 0
double
norm_frob_sq(OperatorComponent const& x)
{
   // TODO: this implementation is not very efficient
   return trace(local_inner_prod(herm(x), x)).real();
   //   return norm_frob_sq(local_inner_prod(herm(x), x));
}
#endif

#if 0
OperatorComponent
shift_left(BasisList const& ThisBasis, BasisList const& IncomingBasis)
{
   // Components Result(a,b)(c,d) = delta_{a,d} delta_{b,c}
   // where a,d are from IncomingBasis, b,c are from ThisBasis
   OperatorComponent Result(IncomingBasis, ThisBasis, ThisBasis, IncomingBasis);

   ProductBasis<BasisList, BasisList>

   SimpleOperator t = SimpleOperator::make_identity(ThisBasis);
   SimpleOperator in = SimpleOperator::make_identity(IncomingBasis);

   QuantumNumber Ident(ThisBasis.GetSymmetryList());
   for (int i = 0; i < IncomingBasis.size(); ++i)
   {
      for (int j = 0; j < ThisBasis.size(); ++j)
      {
         SimpleOperator x(ThisBasis, IncomingBasis, Ident);
         x(j,i) = 1;
      }
   }
}
#endif

OperatorComponent
translate_right(BasisList const& LeftBasis, BasisList const& ThisBasis)
{
   // Components Result(a,b)(c,d)_k = delta_{a,c} delta_{b,d} (2k+1)(2c_1)
   // where a,c are from LeftBasis, b,d are from ThisBasis
   OperatorComponent Result(LeftBasis, ThisBasis, LeftBasis, ThisBasis);

   QuantumNumber Ident(ThisBasis.GetSymmetryList());
   for (unsigned i = 0; i < LeftBasis.size(); ++i)
   {
      for (unsigned j = 0; j < ThisBasis.size(); ++j)
      {
         QuantumNumberList k = inverse_transform_targets(ThisBasis[j], LeftBasis[i]);
         SimpleRedOperator x(LeftBasis, ThisBasis);
         for (auto const& q : k)
         {
            SimpleOperator xComponent(LeftBasis, ThisBasis, q);
            xComponent.insert(i,j, qdim(q) / qdim(ThisBasis[j]));
	    x.insert(std::move(xComponent));
         }
         Result.insert(i,j, x);
      }
   }
   return Result;
}

std::ostream&
operator<<(std::ostream& out, OperatorComponent const& op)
{
   out << "OperatorComponent: Basis1=" << op.Basis1() << "\nBasis2=" << op.Basis2() << '\n';
   for (auto const& r : op)
   {
      for (auto const& c : r)
      {
         out << "element (" << r.row() << "," << c.col() << ") = " << c.value << '\n';
      }
   }
   return out;
}

OperatorComponent
tensor_sum(OperatorComponent const& A, OperatorComponent const& B,
           SumBasis<BasisList> const& B1, SumBasis<BasisList> const& B2)
{
   DEBUG_CHECK_EQUAL(A.LocalBasis1(), B.LocalBasis1());
   DEBUG_CHECK_EQUAL(A.LocalBasis2(), B.LocalBasis2());

   using ::copy;

   OperatorComponent Result(A.LocalBasis1(), A.LocalBasis2(), B1.Basis(), B2.Basis());

   // set the matrix elements corresponding to operator A (subspace 0 in the SumBasis)
   for (auto const& r : A)
   {
      for (auto const& c : r)
      {
         Result.insert(B1(0,r.row()), B2(0,c.col()), copy(c.value));
      }
   }

   // set the matrix elements corresponding to operator A (subspace 1 in the SumBasis)
   for (auto const& r : B)
   {
      for (auto const& c : r)
      {
         Result.insert(B1(1,r.row()), B2(1,c.col()), copy(c.value));
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

   using ::copy;

   OperatorComponent Result(A.LocalBasis1(), A.LocalBasis2(), A.Basis1(), B2.Basis());

   // set the matrix elements corresponding to operator A (subspace 0 in the SumBasis)
   for (auto const& r : A)
   {
      for (auto const& c : r)
      {
         Result.insert(r.row(), B2(0,c.col()), copy(c.value));
      }
   }

   // set the matrix elements corresponding to operator A (subspace 1 in the SumBasis)
   for (auto const& r : B)
   {
      for (auto const& c : r)
      {
         Result.insert(r.row(), B2(1,c.col()), copy(c.value));
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
      return copy(B);
   if (B.is_null())
      return copy(A);
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

   using ::copy;

   OperatorComponent Result(A.LocalBasis1(), A.LocalBasis2(), B1.Basis(), A.Basis2());

   // set the matrix elements corresponding to operator A (subspace 0 in the SumBasis)
   for (auto const& r : A)
   {
      for (auto const& c : r)
      {
         Result.insert(B1(0,r.row()), c.col(), copy(c.value));
      }
   }

   // set the matrix elements corresponding to operator A (subspace 1 in the SumBasis)
   for (auto const& r : B)
   {
      for (auto const& c : r)
      {
         Result.insert(B1(0,r.row()), c.col(), copy(c.value));
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
      return copy(B);
   if (B.is_null())
      return copy(A);
   return tensor_col_sum(A, B, SumBasis<BasisList>(A.Basis1(), B.Basis1()));
}

OperatorComponent prod(OperatorComponent const& A, SimpleOperator const& Op, double Tol)
{
   PRECONDITION(is_scalar(Op.TransformsAs()))(Op.TransformsAs());
   DEBUG_PRECONDITION_EQUAL(A.Basis2(), Op.Basis1());

   OperatorComponent Result(A.LocalBasis1(), A.LocalBasis2(), A.Basis1(), Op.Basis2());
   // for our convention on the coupling constants, this should be a trivial operation
   // of a standard matrix-matrix multiply (the nested operation is a SimpleRedOperator * scalar)
   Result.data() = copy(A.data()) * Op.data();
   return Result;
}

OperatorComponent prod(SimpleOperator const& Op, OperatorComponent const& A, double Tol)
{
   PRECONDITION(is_scalar(Op.TransformsAs()))(Op.TransformsAs());
   DEBUG_PRECONDITION_EQUAL(Op.Basis2(), A.Basis1());

   OperatorComponent Result(A.LocalBasis1(), A.LocalBasis2(), Op.Basis1(), A.Basis2());
   Result.data() = Op.data() * A.data();
   return Result;
}

OperatorComponent prod(OperatorComponent const& A, HermitianProxy<SimpleOperator> const& Op, double Tol)
{
   PRECONDITION(is_scalar(Op.base().TransformsAs()))(Op.base().TransformsAs());
   DEBUG_PRECONDITION_EQUAL(A.Basis2(), Op.base().Basis2());

   OperatorComponent Result(A.LocalBasis1(), A.LocalBasis2(), A.Basis1(), Op.base().Basis1());
   // for our convention on the coupling constants, this should be a trivial operation
   // of a standard matrix-matrix multiply (the nested operation is a SimpleRedOperator * scalar)
   Result.data() = A.data() * herm(Op.base().data());
   return Result;
}

OperatorComponent prod(HermitianProxy<SimpleOperator> const& Op, OperatorComponent const& A, double Tol)
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
            HermitianProxy<SimpleOperator> const& y)
{
   return prod(x, prod(Op, y));
}

OperatorComponent local_tensor_prod(OperatorComponent const& A, OperatorComponent const& B)
{
   DEBUG_PRECONDITION_EQUAL(A.Basis2(), B.Basis1());

   using blas::norm_frob;

   ProductBasis<BasisList, BasisList> LB1(A.LocalBasis1(), B.LocalBasis1());
   ProductBasis<BasisList, BasisList> LB2(A.LocalBasis2(), B.LocalBasis2());

   OperatorComponent Result(LB1.Basis(), LB2.Basis(), A.Basis1(), B.Basis2());

   typedef OperatorComponent::data_type MatType;

   for (auto const& rA : A)
   {
      for (auto const& cA : rA)
      {
	 for (auto const& cB : B.row(cA.col()))
	 {
	    for (auto const& Acomponent : cA.value)
	    {
	       for (auto const& Bcomponent : cB.value)
	       {
		  QuantumNumberList qL = Acomponent.TransformsAs() * Bcomponent.TransformsAs();
		  for (auto const& q : qL)
		  {
                     // Ensure that the final operator is a possible target
                     if (is_transform_target(B.qn2(cB.col()), q, A.qn1(rA.row())))
                     {
                        real Coeff = product_coefficient(Acomponent.TransformsAs(), Bcomponent.TransformsAs(), q,
							 A.qn1(rA.row()), B.qn2(cB.col()), B.qn1(cA.col()));
			// TODO: the product_coefficient should handle small values better and return true zero
                        if (norm_frob(Coeff) > 1E-14)
			   Result.add(rA.row(), cB.col(), Coeff * tensor_prod(Acomponent, Bcomponent, q));
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

// TODO: the old implemention looks buggy
SimpleOperator
local_inner_prod(HermitianProxy<OperatorComponent> const& A, OperatorComponent const& B)
{
   DEBUG_PRECONDITION_EQUAL(A.base().Basis1(), B.Basis1());

   // The coupling coefficient for the auxiliary basis is a simple rescaling,
   // same as for the tensor ScalarProd operation

   SimpleOperator Result(A.base().Basis2(), B.Basis2());

   // AI/BI point at the same quantum number, and is the inner summation of the product
   for (int i = 0; i < B.size1(); ++i)
   {
      double InnerDegree = qdim(B.qn1(i));
      for (auto const& cA : A.base().row(i))
      {
	 for (auto const& cB : B.row(i))
	 {
            if (A.base().qn2(cA.col()) == B.qn2(cB.col()))
               Result.add(cA.col(), cB.col(), InnerDegree / qdim(A.base().qn2(cA.col())) * inner_prod(cA.value, cB.value));
	 }
      }
   }
   Result.debug_check_structure();
   return Result;
}


#if 0
// this doesn't appear to be used
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
#endif

SimpleOperator
local_inner_tensor_prod(HermitianProxy<OperatorComponent> const& A, OperatorComponent const& B)
{
   DEBUG_PRECONDITION_EQUAL(A.base().Basis1(), B.Basis1());

   using blas::norm_frob;

   ProductBasis<BasisList, BasisList> PBasis1(adjoint(A.base().Basis1()), B.Basis1());
   ProductBasis<BasisList, BasisList> PBasis2(adjoint(A.base().Basis2()), B.Basis2());

   QuantumNumbers::QuantumNumber Ident(B.GetSymmetryList());

   SimpleOperator Result(PBasis1.Basis(), PBasis2.Basis(), Ident);

   for (auto const& rA : A.base())
   {
      for (auto const& rB : B)
      {
	 for (auto const& cA : rA)
	 {
	    for (auto const& cB : rB)
	    {
	       for (auto const& b1 : PBasis1(rA.row(), rB.row()))
	       {
		  for (auto const& b2 : PBasis2(cA.col(), cB.col()))
		  {
		     for (auto const& Areduced : cA.value)
		     {
			for (auto const& Breduced : cB.value)
			{
                           // we are taking the inner product, so we require the quantum numbers match
			   if (Areduced.TransformsAs() != Breduced.TransformsAs())
			      continue;

			   if (!is_transform_target(PBasis2[b2], Ident, PBasis1[b1]))
			      continue;

			   auto Coeff = tensor_coefficient(PBasis1, PBasis2,
							   adjoint(Areduced.TransformsAs()), Breduced.TransformsAs(), Ident,
							   rA.row(), rB.row(), b1,
							   cA.col(), cB.col(), b2);
                           if (norm_frob(Coeff) > 1E-14)
                           {
			      Result.add(b1, b2, Coeff * inner_prod(Areduced, Breduced));
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
RealSimpleOperator
local_norm_frob_sq(OperatorComponent const& A)
{
   RealSimpleOperator Result(A.Basis1(), A.Basis2());
   for (OperatorComponent::const_iterator i = iterate(A); i; ++i)
   {
      for (OperatorComponent::const_inner_iterator j = iterate(i); j; ++j)
      {
         Result(i,j) = norm_frob_sq(*j);
      }
   }
   return Result;
}
#endif

#if 0
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
#endif

OperatorComponent aux_tensor_prod(OperatorComponent const& ML, OperatorComponent const& MR)
{
   DEBUG_PRECONDITION_EQUAL(ML.LocalBasis2(), MR.LocalBasis1());
   ProductBasis<BasisList, BasisList> PBasis1(ML.Basis1(), MR.Basis1());
   ProductBasis<BasisList, BasisList> PBasis2(ML.Basis2(), MR.Basis2());

   using blas::norm_frob;

   OperatorComponent Result(ML.LocalBasis1(), MR.LocalBasis2(), PBasis1.Basis(), PBasis2.Basis());

   for (auto const& rML : ML)
   {
      for (auto const& cML : rML)
      {
	 for (auto const& rMR : MR)
	 {
	    for (auto  const& cMR : rMR)
	    {

	       for (auto const& b1 : PBasis1(rML.row(), rMR.row()))
	       {
		  for (auto const& b2 : PBasis2(cML.col(), cMR.col()))
		  {
		     for (auto const& Areduced : cML.value)
		     {
			for (auto const& Breduced : cMR.value)
			{
			   for (auto const& q : transform_targets(Areduced.TransformsAs(), Breduced.TransformsAs()))
			   {
			      if (is_transform_target(PBasis2[b2], q, PBasis2[b1]))
			      {
                                 double Coeff = tensor_coefficient(PBasis1, PBasis2,
                                                                   Areduced.TransformsAs(), Breduced.TransformsAs(), q,
                                                                   rML.row(), rMR.row(), b1,
                                                                   cML.col(), cMR.col(), b2);
                                 if (norm_frob(Coeff) > 1E-14)
                                 {
				    Result.add(b1, b2, Coeff * prod(Areduced, Breduced, q));
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
// Result[s'][i'a,j'b] = Op(s',s)(a,b) * S(s)(i,j)
StateComponent aux_tensor_prod(OperatorComponent const& Op, StateComponent const& S)
{
   DEBUG_PRECONDITION_EQUAL(Op.LocalBasis2(), S.LocalBasis());
   ProductBasis<BasisList, VectorBasis> PBasis1(Op.Basis1(), S.Basis1());
   ProductBasis<BasisList, VectorBasis> PBasis2(Op.Basis2(), S.Basis2());

   StateComponent Result(Op.LocalBasis1(), PBasis1.Basis(), PBasis2.Basis());

   // iterate over the aux indices of the MPO
   for (OperatorComponent::const_iterator iOp = iterate(Op); iOp; ++iOp)
   {
      for (OperatorComponent::const_inner_iterator jOp = iterate(iOp); jOp; ++jOp)
      {
         // iterate over the reducible components at MPO(i,j)
         for (SimpleRedOperator::const_iterator RedOp = jOp->begin(); RedOp != jOp->end(); ++RedOp)
         {
            // iterate over the matrix elements, s' is jInner.index1(), s is jInner.index2()
            for (SimpleOperator::const_iterator iInner = iterate(*RedOp); iInner; ++iInner)
            {
               for (SimpleOperator::const_inner_iterator jInner = iterate(iInner); jInner; ++jInner)
               {
                  // Now iterate over the vector basis of the MPS
                  for (MatrixOperator::const_iterator iS = iterate(S[jInner.index2()]); iS; ++iS)
                  {
                     for (MatrixOperator::const_inner_iterator jS = iterate(iS); jS; ++jS)
                     {
                        // now we have all the indices, map them onto the final index via the
                        // product basis

                        ProductBasis<BasisList, VectorBasis>::const_iterator TiIter,
                           TiEnd = PBasis1.end(jOp.index1(), jS.index1());
                        for (TiIter = PBasis1.begin(jOp.index1(), jS.index1()); TiIter != TiEnd; ++TiIter)
                        {
                           ProductBasis<BasisList, VectorBasis>::const_iterator TjIter,
                              TjEnd = PBasis2.end(jOp.index2(), jS.index2());
                           for (TjIter = PBasis2.begin(jOp.index2(), jS.index2()); TjIter != TjEnd;
                                ++TjIter)
                           {
                              if (!is_transform_target(PBasis2[*TjIter], Result.LocalBasis()[jInner.index1()],
                                                      PBasis1[*TiIter]))
                                 continue;

                              double Coeff = tensor_coefficient(PBasis1, PBasis2,
                                                                RedOp->TransformsAs(),
                                                                S.LocalBasis()[jInner.index2()],
                                                                Result.LocalBasis()[jInner.index1()],
                                                                jOp.index1(), jS.index1(), *TiIter,
                                                                jOp.index2(), jS.index2(), *TjIter);
                              if (LinearAlgebra::norm_2(Coeff) > 1E-14)
                              {
                                 add_element(Result[jInner.index1()].data(), *TiIter, *TjIter,
                                             conj_phase(S.LocalBasis()[jInner.index2()],
                                                        RedOp->TransformsAs(),
                                                        Result.LocalBasis()[jInner.index1()])
                                             * Coeff * (*jInner) * (*jS)
                                             );
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
#endif

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

#if 0
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
#endif

SimpleOperator TruncateBasis1(OperatorComponent& A)
{
   // We want to work from the last row to the first, so that we preserve the last row exactly.
   // We don't have to worry about accidentally eliminating the first row of a triangular MPO,
   // since if it is triangular then the last column cannot be parallel to anything else

   // iterate through the rows of A and remove any that are parallel with a previous.
   // We cannot safely remove rows that are merely linearly dependent on the previous
   // rows without potential trouble with catastrophic cancellation.

   using blas::norm_frob;

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
      if (imat <= Scale * RemoveRowEpsilon)
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
         if (norm_frob(ijmat) / std::sqrt(imat*jmat) > (1.0 - TruncateOverlapEpsilon))
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

   A = std::move(tA);
   return Reg;
}

SimpleOperator TruncateBasis2(OperatorComponent& A)
{
   // We want to work from the first column to last, so that we preserve the first column exactly.
   // For a BasicTriangularMPO the first column will contain the identity.
   // We don't have to worry about accidentally eliminating the last column of a triangular MPO,
   // since if it is triangular then the last column cannot be parallel to anything else.

   using blas::norm_frob;

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
      if (imat <= Scale * RemoveRowEpsilon)
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
         if (norm_frob(ijmat) / std::sqrt(imat*jmat) > (1.0 - TruncateOverlapEpsilon))
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

   A = std::move(tA);
   return Reg;
}

#if 0
SimpleOperator TruncateBasis1MkII(OperatorComponent& A, double Epsilon)
{
   // We want to work from the last row to the first, so that we preserve the last row exacrly.
   // For a BasicTriangularMPO the last row will contain the identity.
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
   double Scale = std::sqrt(norm_frob_sq(A) / (A.Basis1().total_degree() * A.Basis2().total_degree()));

   //   TRACE(Scale)(norm_frob_sq(A))(A.Basis1().total_degree())(A.Basis2().total_degree());
   //   TRACE(A);
   TRACE(Scale)(norm_frob_sq(A))(A.Basis1().total_degree())(A.Basis2().total_degree());
   TRACE(A);

   // make a dense matrix
   blas::Matrix<SimpleRedOperator> M(A.data());

   // these are stored in the reverse order
   std::vector<blas::Vector<SimpleRedOperator> > Rows;
   std::vector<blas::Vector<std::complex<double> > > T;

   std::vector<double> RowNormSq;

   std::vector<QuantumNumbers::QuantumNumber> NewBasis1Q;

   //   double Normalization = A.LocalBasis2().total_degree();

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
         //      Rows.push_back((std::sqrt(Normalization / NextRowNormSq)) * NextRow);
         Rows.push_back(NextRow);
         T.push_back(LinearAlgebra::Vector<std::complex<double> >(M.size1(), 0.0));
         //T[Rows.size()-1][r] = std::sqrt(NextRowNormSq / Normalization);
         T[Rows.size()-1][r] = 1;
         NewBasis1Q.push_back(A.Basis1()[r]);
         RowNormSq.push_back(NextRowNormSq);
         //RowNormSq.push_back(Normalization);  // since we've already normalized it
      }
      else
      {
         TRACE(NextRowNormSq)(Epsilon);
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
   CHECK(norm_frob(A - ACheck) <= Scale*QREpsilon)(A-ACheck)(Scale)(TruncateOverlapEpsilon);

   ANew.check_structure();
   Trunc.check_structure();

   A = ANew;
   return Trunc;
}
#endif

#if 0
SimpleOperator TruncateBasis2MkII(OperatorComponent& A, double Epsilon)
{
   // We want to work from the first column to last, so that we preserve the first column exactly.
   // For a BasicTriangularMPO the first column will contain the identity.
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
   double Scale = std::sqrt(norm_frob_sq(A) / (A.Basis1().total_degree() * A.Basis2().total_degree()));

   // make a dense matrix
   LinearAlgebra::Matrix<SimpleRedOperator> M = A.data();

   std::vector<LinearAlgebra::Vector<SimpleRedOperator> > Columns;
   std::vector<LinearAlgebra::Vector<std::complex<double> > > T;

   std::vector<double> ColNormSq;

   BasisList NewBasis2(A.Basis2().GetSymmetryList());

   //double Normalization = A.LocalBasis1().total_degree();

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
         //Columns.push_back((std::sqrt(Normalization / NextColNormSq)) * NextCol);
         Columns.push_back(NextCol);
         T.push_back(LinearAlgebra::Vector<std::complex<double> >(M.size2(), 0.0));
         //T[Columns.size()-1][c] = std::sqrt(NextColNormSq / Normalization);
         T[Columns.size()-1][c] = 1;
         NewBasis2.push_back(A.Basis2()[c]);
         //ColNormSq.push_back(Normalization);  // since we've already normalized it
         ColNormSq.push_back(NextColNormSq);
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
   CHECK(norm_frob(A - ACheck) <= Scale*QREpsilon)(norm_frob(A - ACheck))(Scale);

   ANew.check_structure();
   Trunc.check_structure();

   A = ANew;
   return Trunc;
}
#endif

// compress the bond dimension of basis 2 using the stabilized
// 'generalized deparallelization' algorithm

#if 0
SimpleOperator CompressBasis2_LinDep(OperatorComponent& A)
{
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

   // Get the squared norm of each component of the operator
   RealSimpleOperator const LocalNorms = local_norm_frob_sq(A);

   // The transform matrix.  Trans[c] is the transform matrix for column c.
   // If Trans[c] is empty, then the corresponding column is linearly dependent and will be removed.
   // Effecively, Trans[i][j] is the entry (i,j) of the transform matrix of A = A' * Trans
   std::vector<LinearAlgebra::MapVector<std::complex<double>>> Trans(A.Basis2().size());

   for (int i = 0; i < int(A.Basis2().size()); ++i)
   {
      // Is the column empty?
      bool Empty = true;
      for (int r = 0; r < int(A.Basis1().size()); ++r)
      {
         if (LocalNorms(r,i) != 0.0)
         {
            Empty = false;
            break;
         }
      }
      if (Empty)
      {
         // We don't need to do anything here, since Trans[i] is empty that signals that the
         // column is deleted.
         continue;
      }

      //bool Dependent = false;  // is vector i linearly dependent on some other set?
      // assemble the list of candidate columns that could give a linear dependence
      std::vector<int> Candidates;
      for (int j = 0; j < i; ++j)
      {
         // to be a candidate, every row of column j that is non-zero
         // must also be non-zero in column i.  It must also have the same quantum number.
         if (Trans[j].is_zero())
            continue;

         if (A.Basis2()[j] != A.Basis2()[i])
            continue;

         bool Candidate = true;
         for (int r = 0; i < int(A.Basis1().size()); ++r)
         {
            if (LocalNorms(r,j) > 0 && LocalNorms(r,i) == 0)
            {
               Candidate = false;
               break;
            }
         }
         if (Candidate)
            Candidates.push_back(j);
      }

      // Do we have a possible linear dependency?
      if (!Candidates.empty())
      {
         // Form the matrix of overlaps
         // Need to linearize the local operators into vectors, effectively 'reshaping' the matrix
         std::vector<int> UsedRows;
         std::vector<int> OffsetOfRow;
         int TotalRows = 0;
         for (int r = 0; r < int(A.Basis1().size()); ++r)
         {
            if (LocalNorms(r,i) != 0.0)
            {
               UsedRows.push_back(r);
               OffsetOfRow.push_back(TotalRows);
               int Rows = linear_dimension(A(r,i),
                                           transform_targets(A.Basis2()[i],
                                                            adjoint(A.Basis1()[r])));
               TotalRows += Rows;
            }
         }
         OffsetOfRow.push_back(TotalRows);

         LinearAlgebra::Matrix<std::complex<double>> X(TotalRows, Candidates.size(), 0.0);
         LinearAlgebra::Vector<std::complex<double>> RHS(TotalRows, 0.0);
         for (int r = 0; r < int(UsedRows.size()); ++r)
         {
            double Scale = LocalNorms(r,i);
            for (int c = 0; c < int(Candidates.size()); ++c)
            {
               X(LinearAlgebra::range(OffsetOfRow[r], OffsetOfRow[r+1]),LinearAlgebra::all)(LinearAlgebra::all, c)
                  = (1.0/Scale) * linearize(A(UsedRows[r],Candidates[c]),
                                            transform_targets(A.Basis2()[Candidates[c]],
                                                              adjoint(A.Basis1()[UsedRows[r]])));
            }
            RHS[LinearAlgebra::range(OffsetOfRow[r], OffsetOfRow[r+1])]
               = (1.0/Scale) * linearize(A(UsedRows[r],i),
                                         transform_targets(A.Basis2()[i], adjoint(A.Basis1()[r])));
         }

         LinearAlgebra::Vector<std::complex<double>> x;
         double Resid;
         std::tie(Resid, x) = LeastSquaresRegularized(X, RHS);

         if (Resid < 1E-14)
         {
            // linear dependency
            for (int c = 0; c < int(size(x)); ++c)
            {
               Trans[Candidates[c]][i] = x[c];
            }
            // next row
            continue;
         }
      }
      else
      {
         // the column is linearly independent
         Trans[i] = LinearAlgebra::MapVector<std::complex<double>>(A.Basis2().size());
         Trans[i][i] = 1.0;
      }
   }


#if 0
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
#endif
}
#endif

#if 0
// See optimized version in f-optim.cpp
StateComponent
operator_prod(OperatorComponent const& M,
              StateComponent const& A,
              StateComponent const& F,
              HermitianProxy<StateComponent> const& B)
{
   PRECONDITION_EQUAL(M.LocalBasis1(), A.LocalBasis());
   PRECONDITION_EQUAL(M.LocalBasis2(), B.base().LocalBasis());
   DEBUG_PRECONDITION_EQUAL(M.Basis2(), F.LocalBasis());
   DEBUG_PRECONDITION_EQUAL(A.Basis2(), F.Basis1());
   DEBUG_PRECONDITION_EQUAL(F.Basis2(), B.base().Basis2());

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
                                  F[J.index2()],
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
#endif

#if 0
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
#endif

#if defined(OLD_OPERATOR_PROD)
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
#endif

#if 0
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

// Result[s'](i',i) = M(s',s)[a',a] herm(E[a'](i',j')) herm(B[s](j',j)) F[a](i,j)
StateComponent
operator_prod_inner(OperatorComponent const& M,
                    HermitianProxy<StateComponent> const& E,
                    HermitianProxy<StateComponent> const& B,
                    StateComponent const& F)
{
   DEBUG_PRECONDITION_EQUAL(M.Basis1(), E.base().LocalBasis());
   DEBUG_PRECONDITION_EQUAL(M.Basis2(), F.LocalBasis());
   DEBUG_PRECONDITION_EQUAL(M.LocalBasis2(), A.LocalBasis());

   DEBUG_PRECONDITION_EQUAL(A.Basis2(), F.Basis1());
   DEBUG_PRECONDITION_EQUAL(E.base().Basis1(), A.Basis1());

   StateComponent Result(M.LocalBasis1(), E.base().Basis2(), F.Basis2());

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

                  // iterate over


                  Result[S.index1()] += (*S) * triple_prod(herm(E.base()[J.index1()]),
                                                           A[S.index2()],
                                                           F[J.index2()],
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
#endif

OperatorComponent
project_rows(OperatorComponent const& x, std::set<int> const& Rows)
{
   using ::copy;
   
   BasisList ProjectedBasis(x.Basis1().GetSymmetryList());
   std::vector<int> Map(x.Basis1().size(), -1);
   for (std::set<int>::const_iterator I = Rows.begin(); I != Rows.end(); ++I)
   {
      Map[*I] = ProjectedBasis.size();
      ProjectedBasis.push_back(x.Basis1()[*I]);
   }

   OperatorComponent Result(x.LocalBasis1(), x.LocalBasis2(), ProjectedBasis, x.Basis2());
   for (auto const& r : x)
   {
      if (Map[r.row()] < 0)
	 continue;

      for (auto const& c : r)
      {
	 Result.insert(Map[r.row()], c.col(), copy(c.value));
      }
   }

   return Result;
}

// project onto the given columns of the component
OperatorComponent
project_columns(OperatorComponent const& x, std::set<int> const& Cols)
{
   using ::copy;

   BasisList ProjectedBasis(x.Basis2().GetSymmetryList());
   std::vector<int> Map(x.Basis2().size(), -1);
   for (std::set<int>::const_iterator I = Cols.begin(); I != Cols.end(); ++I)
   {
      Map[*I] = ProjectedBasis.size();
      ProjectedBasis.push_back(x.Basis2()[*I]);
   }

   OperatorComponent Result(x.LocalBasis1(), x.LocalBasis2(), x.Basis1(), ProjectedBasis);
   for (auto const& r : x)
   {
      for (auto const& c : r)
      {
	 if (Map[c.col()] >= 0)
	 {
	    Result.insert(r.row(), Map[c.col()], copy(c.value));
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
   typedef std::tuple<QuantumNumbers::QuantumNumber, int, int> PartialProdIndex;
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
   blas::Matrix<std::complex<double> > Mat(LeftIndex.size(), RightIndex.size(), 0.0);
   for (PartialProdType::const_iterator I = PartialProd.begin(); I != PartialProd.end(); ++I)
   {
      PartialProdIndex Left(I->first.qLeft, I->first.Left1, I->first.Left2);
      PartialProdIndex Right(I->first.qRight, I->first.Right1, I->first.Right2);
      Mat(LeftIndex[Left], RightIndex[Right]) = I->second;
   }

   // Perform the singular decomposition
   int mn = std::min(Mat.rows(), Mat.cols());
   blas::Matrix<complex> U(Mat.rows(), mn), Vh(mn, Mat.cols());
   blas::Vector<double> D(mn);
   SingularValueDecomposition(std::move(Mat), U, D, Vh);

   // This sets the scale of the singular values.  Any singular values smaller than
   // KeepThreshold are removed.
   double KeepThreshold = std::numeric_limits<double>::epsilon() * sum(D) * 10;

   // split into pairs of operators, keeping only the non-zero singular values
   std::vector<std::pair<SimpleRedOperator, SimpleRedOperator> > Result;
   for (unsigned k = 0; k < D.size(); ++k)
   {
      DEBUG_TRACE(D[k]);
      if (D[k] > KeepThreshold)
      {
         double Coeff = std::sqrt(D[k]);  // distribute sqrt(D) to each operator
         SimpleRedOperator L(B1.Left(), B2.Left());
         for (unsigned x = 0; x < U.rows(); ++x)
         {
            if (norm_frob(U(x,k)) > std::numeric_limits<double>::epsilon() * 10)
            {
               L.project(std::get<0>(LeftReverse[x]))
                  (std::get<1>(LeftReverse[x]), std::get<2>(LeftReverse[x])) = U(x,k) * Coeff;
            }
         }

         SimpleRedOperator R(B1.Right(), B2.Right());
         for (unsigned x = 0; x < Vh.cols(); ++x)
         {
            if (norm_frob(Vh(k,x)) > std::numeric_limits<double>::epsilon() * 10)
            {
               R.project(std::get<0>(RightReverse[x]))
                  (std::get<1>(RightReverse[x]), std::get<2>(RightReverse[x])) = Vh(k,x) * Coeff;
            }
         }
         Result.emplace_back(std::move(L),std::move(R));
      }
   }

#if !defined(NDEBUG)
   SimpleOperator OpTest;
   for (unsigned i = 0; i < Result.size(); ++i)
   {
      // We really need an IrredProd for reducible tensors, it would be more efficient here
      OpTest += project(Result[i].first * Result[i].second, Op.TransformsAs());
   }
   CHECK(norm_frob(OpTest - Op) < 1E-13 * norm_frob(OpTest))(Op)(OpTest)(norm_frob(OpTest-Op));
#endif

   return Result;
}

std::pair<OperatorComponent, OperatorComponent>
decompose_local_tensor_prod(OperatorComponent const& Op,
                            ProductBasis<BasisList, BasisList> const& B1,
                            ProductBasis<BasisList, BasisList> const& B2)
{
   // A scale factor for determining whether to keep numerically small values
   //double ScaleFactor = norm_frob_sq(Op) / std::sqrt(Op.Basis1().size() * Op.Basis2().size() * Op.LocalBasis1().size() * Op.LocalBasis2().size());

   typedef std::map<Tensor::PartialProdIndex, std::complex<double> > PartialProdType;

   // ComponentIndex is an index to a particular matrix element of the MPO transforming as q,
   // between local indices local1,local2 and auxiliary index aux
   // ComponentIndex(aux, local1, local2, q)
   typedef std::tuple<int, int, int, QuantumNumbers::QuantumNumber> ComponentIndex;

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

   for (auto const& rOp : Op)
   {
      int LeftAux1 = rOp.row();
      for (auto const& cOp : rOp)
      {
	 int RightAux2 = cOp.col();
	 for (auto const& Comp : cOp.value)
	 {
            // Decompose the tensor product component
            PartialProdType PartialProd = Tensor::decompose_tensor_prod(Comp, B1, B2);

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
	    for (auto const& P : PartialProd)
	    {
	       for (auto const& Qpp :  Op.qn2(cOp.col()) * P.first.qRight)
	       {
		  // Skip quantum numbers that are not possible (structural zero in the coupling coefficient)
                  if (!is_transform_target(Qpp,  P.first.qLeft, Op.qn1(rOp.row())))
		     continue;
		      
		  // Qpp is the quantum number of the inner basis
		  
		  auto Coeff = inverse_product_coefficient(P.first.qLeft, P.first.qRight, Comp.TransformsAs(),
                                                             Op.qn1(rOp.row()), Op.qn2(cOp.col()), Qpp);
                  if (norm_frob(Coeff) > 1E-14)
                  {
                     ComponentIndex LeftIndex
                        = std::make_tuple(LeftAux1, P.first.Left1, P.first.Left2, P.first.qLeft);
                     if (LeftMapping[Qpp].find(LeftIndex) == LeftMapping[Qpp].end())
                     {
                        LeftMapping[Qpp][LeftIndex] = LeftRevMapping[Qpp].size();
                        LeftRevMapping[Qpp].push_back(LeftIndex);
                     }

                     ComponentIndex RightIndex
                        = std::make_tuple(RightAux2, P.first.Right1, P.first.Right2, P.first.qRight);
                     if (RightMapping[Qpp].find(RightIndex) == RightMapping[Qpp].end())
                     {
                        RightMapping[Qpp][RightIndex] = RightRevMapping[Qpp].size();
                        RightRevMapping[Qpp].push_back(RightIndex);
                     }

                     DecomposedOperator[Qpp][std::make_pair(LeftIndex, RightIndex)]
                        = Coeff * P.second;
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
      blas::Matrix<complex> Mat(LeftMapping[Qpp].size(), RightMapping[Qpp].size(), 0.0);
      for (DecomposedMPOType::const_iterator I = Q->second.begin(); I != Q->second.end(); ++I)
      {
         Mat(LeftMapping[Qpp][I->first.first], RightMapping[Qpp][I->first.second]) = I->second;
      }

      // Do the SVD
      int mn = std::min(Mat.rows(), Mat.cols());
      blas::Matrix<std::complex<double> > U(Mat.rows(), mn), Vh(mn, Mat.cols());
      blas::Vector<double> D(mn);
      SingularValueDecomposition(std::move(Mat), U, D, Vh);

      DEBUG_TRACE(D);

      // TODO: this used to be 'max'
      double m = sum(D);
      if (m > LargestSingular)
         LargestSingular = m;
   }

   // This sets the scale of the singular values.  Any singular values smaller than
   // KeepThreshold are removed.  The coefficient here is determined as appropriate
   // for culling unwanted matrix elements from exponentials.  Possibly it needs to scale
   // with be basis size or the number of singular values in some way?
   double KeepThreshold = std::numeric_limits<double>::epsilon() * LargestSingular * 20;

   DEBUG_TRACE(KeepThreshold);

   // Second pass, collate the kept singular values
   for (DecomposedByQuantumType::const_iterator Q = DecomposedOperator.begin();
        Q != DecomposedOperator.end(); ++Q)
   {
      QuantumNumbers::QuantumNumber Qpp = Q->first;

      // assemble the matrix
      blas::Matrix<std::complex<double> > Mat(LeftMapping[Qpp].size(), RightMapping[Qpp].size(), 0.0);
      for (DecomposedMPOType::const_iterator I = Q->second.begin(); I != Q->second.end(); ++I)
      {
         Mat(LeftMapping[Qpp][I->first.first], RightMapping[Qpp][I->first.second]) = I->second;
      }

      // Do the SVD
      int mn = std::min(Mat.rows(), Mat.cols());
      blas::Matrix<std::complex<double> > U(Mat.rows(), mn), Vh(mn, Mat.cols());
      blas::Vector<double> D(mn);
      SingularValueDecomposition(std::move(Mat), U, D, Vh);

      // Assemble the singular vectors
      for (unsigned k = 0; k < D.size(); ++k)
      {
         // Skip over vectors that don't reach the threshold
         if (D[k] < KeepThreshold)
            continue;

         // Distribute the singular value to the two sides.  It doesn't matter much how we do this;
         // sqrt(D) to each side would work fine.  But we get better balance by distributing according
         // to the dimenions of the left/right local Hilbert spaces
         double LeftCoeff = std::pow(D[k], log(B1.Left().total_degree()) / (log(B1.Left().total_degree()) + log(B1.Right().total_degree())));
         double RightCoeff = D[k] / LeftCoeff;

         SingularVectorType LeftVec;
         for (unsigned x = 0; x < U.rows(); ++x)
         {
            // We might have some components that are in forbidden quantum number sectors - these should be small.
            // In the auxiliary space, we should have a matrix element that transforms as q,
            // with q1 = aux, q2 = Qpp
            // In the local space, we have a matrix element that transforms as q,
            // with q1 = local1, q2 = local2,
            // where
            // LeftRevMapping[Qpp][x] is a ComponentIndex(aux, local1, local2, q)

            if (!is_transform_target(Qpp, std::get<3>(LeftRevMapping[Qpp][x]),
                                     Op.Basis1()[std::get<0>(LeftRevMapping[Qpp][x])])
                || !is_transform_target(B2.Left()[std::get<2>(LeftRevMapping[Qpp][x])],
                                        std::get<3>(LeftRevMapping[Qpp][x]),
                                        B1.Left()[std::get<1>(LeftRevMapping[Qpp][x])]))
            {
               TRACE("Ignoring forbidden off-diagonal matrix element")
                  (U(x,k))(Qpp)(std::get<3>(LeftRevMapping[Qpp][x]))
                  (Op.Basis1()[std::get<0>(LeftRevMapping[Qpp][x])])
                  (B2.Left()[std::get<2>(LeftRevMapping[Qpp][x])])
                  (Qpp)(B1.Left()[std::get<1>(LeftRevMapping[Qpp][x])]);
            }
            else if (norm_frob(U(x,k)) > std::numeric_limits<double>::epsilon() * 10)
            {
               LeftVec[LeftRevMapping[Qpp][x]] = LeftCoeff*U(x,k);
            }
         }

         SingularVectorType RightVec;
         for (unsigned x = 0; x < Vh.cols(); ++x)
         {
            if (norm_frob(Vh(k,x)) > std::numeric_limits<double>::epsilon() * 10)
            {
               RightVec[RightRevMapping[Qpp][x]] = RightCoeff*Vh(k,x);
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

   //   static_assert(std::is_nothrow_move_constructible<OperatorComponent>::value, "");

   return std::pair<OperatorComponent,OperatorComponent>(std::move(MA), std::move(MB));

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
	    SimpleOperator x(B1.Left(), B2.Left(), std::get<3>(I->first));
	    x.insert(std::get<1>(I->first), std::get<2>(I->first), I->second);
	    
	    MA.add(std::get<0>(I->first), b, x);
	 }

         SingularVectorType RightVec = RightSingularVectors[Qpp][i];
         for (SingularVectorType::const_iterator I = RightVec.begin(); I != RightVec.end(); ++I)
         {
	    SimpleOperator x(B1.Right(), B2.Right(), std::get<3>(I->first));
	    x.insert(std::get<1>(I->first), std::get<2>(I->first), I->second);
	    
	    MB.add(b, std::get<0>(I->first), x);
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

   return std::make_pair(std::move(MA), std::move(MB));
}

std::pair<OperatorComponent, OperatorComponent>
decompose_local_tensor_prod(SimpleOperator const& Op,
                            ProductBasis<BasisList, BasisList> const& B1,
                            ProductBasis<BasisList, BasisList> const& B2)
{
   OperatorComponent OpC(Op.Basis1(), Op.Basis2(),
                         make_single_basis(Op.TransformsAs()), make_vacuum_basis(Op.GetSymmetryList()));
   OpC.insert(0, 0, copy(Op));
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

#if 0
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
#endif

void
inplace_conj(OperatorComponent& x)
{
   for (auto& r : x)
   {
      for (auto&& c : r)
      {
	 inplace_conj(c.value);
      }
   }
}

OperatorComponent
conj(OperatorComponent x)
{
   inplace_conj(x);
   return x;
}

OperatorComponent
flip_conj(OperatorComponent const& x)
{
   OperatorComponent Result(adjoint(x.LocalBasis1()), adjoint(x.LocalBasis2()), adjoint(x.Basis1()), adjoint(x.Basis2()));
   for (auto& r : Result)
   {
      for (auto&& c : r)
      {
	 c.value = flip_conj(c.value);
      }
   }
   return Result;
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

void
update_mask_basis1(std::vector<int>& Mask1, SimpleOperator const& Op, std::vector<int> const& Mask2)
{
   for (auto const& r : Op)
   {
      for (auto const& c : r)
      {
	 if (Mask2[c.col()])
	    Mask1[r.row()] = true;
      }
   }
}

void
update_mask_basis1(std::vector<int>& Mask1, SimpleRedOperator const& Op, std::vector<int> const& Mask2)
{
   for (auto const& c : Op)
   {
      update_mask_basis1(Mask1, c, Mask2);
   }
}

void
update_mask_basis1(std::vector<int>& Mask1, OperatorComponent const& Op, std::vector<int> const& Mask2)
{
   for (auto const& r : Op)
   {
      for (auto const& c : r)
      {
         update_mask_basis1(Mask1, c.value, Mask2);
      }
   }
}

void
update_mask_basis2(std::vector<int> const& Mask1, SimpleOperator const& Op, std::vector<int>& Mask2)
{
   for (auto const& r : Op)
   {
      for (auto const& c : r)
      {
	 if (Mask1[r.row()])
	    Mask2[c.col()] = true;
      }
   }
}

void
update_mask_basis2(std::vector<int> const& Mask1, SimpleRedOperator const& Op, std::vector<int>& Mask2)
{
   for (auto const& c : Op)
   {
      update_mask_basis2(Mask1, c, Mask2);
   }
}

void
update_mask_basis2(std::vector<int> const& Mask1, OperatorComponent const& Op, std::vector<int>& Mask2)
{
   for (auto const& r : Op)
   {
      for (auto const& c : r)
      {
         update_mask_basis2(Mask1, c.value, Mask2);
      }
   }
}
