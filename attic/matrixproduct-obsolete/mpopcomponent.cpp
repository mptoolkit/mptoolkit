// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// matrixproduct-obsolete/mpopcomponent.cpp
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

#include "mpopcomponent.h"
#include <tuple>
#include "density.h"
#include "linearalgebra/eigen.h"
#include "common/environment.h"
#include <map>

// the epsilon^2 at which we want to remove matrix elements from operator matrices
// double const OperatorTruncationEpsilon = 0; // for old behaviour
double const OperatorTruncationEpsilon = getenv_or_default("MP_OPERATOREPSILON",
                                                           std::numeric_limits<double>::epsilon());

template <>
void MPOpComponent::check_structure() const
{
   for (MPOpComponent::const_iterator I = this->begin(); I != this->end(); ++I)
   {
      CHECK_EQUAL(I->first, I->second.TransformsAs());
      CHECK_EQUAL(I->second.Basis1(), this->SiteBasis());
      CHECK_EQUAL(I->second.Basis2(), this->SiteBasis());

      I->second.check_structure();

      for (LinearAlgebra::const_iterator<mapped_type>::type X = iterate(I->second); X; ++X)
      {
         for (LinearAlgebra::const_inner_iterator<mapped_type>::type Y = iterate(X); Y; ++Y)
         {
            CHECK_EQUAL(I->first, Y->TransformsAs());
            CHECK_EQUAL(Y->Basis1(), this->Basis1());
            CHECK_EQUAL(Y->Basis2(), this->Basis2());

            Y->check_structure();
         }
      }
   }
}

// helper functor to return the tensor_sum(x, Op), where
// Op is fixed.
struct ApplyTensorSumBind2nd
{
   typedef SimpleOperator const& argument_type;
   typedef SimpleOperator result_type;

   ApplyTensorSumBind2nd(SumBasis<BasisList> const& B1, SumBasis<BasisList> const& B2,
                         SimpleOperator const& Op) : B1_(B1), B2_(B2), Op_(Op) {}

   result_type operator()(argument_type x) const
   {
      return tensor_sum(x, Op_, B1_, B2_);
   }

   SumBasis<BasisList> B1_;
   SumBasis<BasisList> B2_;
   SimpleOperator Op_;
};

struct ApplyTensorSumBind1st
{
   typedef SimpleOperator const& argument_type;
   typedef SimpleOperator result_type;

   ApplyTensorSumBind1st(SumBasis<BasisList> const& B1, SumBasis<BasisList> const& B2,
                         SimpleOperator const& Op) : B1_(B1), B2_(B2), Op_(Op) {}

   result_type operator()(argument_type x) const
   {
      return tensor_sum(Op_, x, B1_, B2_);
   }

   SumBasis<BasisList> B1_;
   SumBasis<BasisList> B2_;
   SimpleOperator Op_;
};

MPOpComponent tensor_sum(MPOpComponent const& A, MPOpComponent const& B,
                         SumBasis<BasisList> const& B1, SumBasis<BasisList> const& B2)
{
   PRECONDITION_EQUAL(A.SiteBasis(), B.SiteBasis());

   typedef MPOpComponent::mapped_type OperatorType;

   MPOpComponent Result(A.SiteBasis(), B1, B2);

   // loop over components in A and add A \oplus B_zero.
   // then loop over components in B and add A_zero \oplus B
   for (MPOpComponent::const_iterator I = A.begin(); I != A.end(); ++I)
   {
      // make a zero operator over the B basis that transforms as the current component
      MPOpComponent::OperatorType BZero(B.Basis1(), B.Basis2(), I->first);
      Result[I->first] = transform(I->second, ApplyTensorSumBind2nd(B1, B2, BZero));
   }

   // do the same for B
   for (MPOpComponent::const_iterator I = B.begin(); I != B.end(); ++I)
   {
      // make a zero operator over the A basis that transforms as the current component
      MPOpComponent::OperatorType AZero(A.Basis1(), A.Basis2(), I->first);
      Result[I->first] += transform(I->second, ApplyTensorSumBind1st(B1, B2, AZero));
   }

   Result.debug_check_structure();
   return Result;
}

struct ApplyTensorRowSumBind2nd
{
   typedef SimpleOperator const& argument_type;
   typedef SimpleOperator result_type;

   ApplyTensorRowSumBind2nd(SumBasis<BasisList> const& B2,
                            SimpleOperator const& Op) : B2_(B2), Op_(Op) {}

   result_type operator()(argument_type x) const
   {
      return tensor_row_sum(x, Op_, B2_);
   }

   SumBasis<BasisList> B2_;
   SimpleOperator Op_;
};

struct ApplyTensorRowSumBind1st
{
   typedef SimpleOperator const& argument_type;
   typedef SimpleOperator result_type;

   ApplyTensorRowSumBind1st(SumBasis<BasisList> const& B2,
                            SimpleOperator const& Op) : B2_(B2), Op_(Op) {}

   result_type operator()(argument_type x) const
   {
      return tensor_row_sum(Op_, x, B2_);
   }

   SumBasis<BasisList> B2_;
   SimpleOperator Op_;
};

MPOpComponent tensor_row_sum(MPOpComponent const& A,
                             MPOpComponent const& B,
                             SumBasis<BasisList> const& B2)
{
   DEBUG_PRECONDITION_EQUAL(A.SiteBasis(), B.SiteBasis());
   DEBUG_PRECONDITION_EQUAL(A.Basis1(), B.Basis1());

   typedef MPOpComponent::mapped_type OperatorType;

   MPOpComponent Result(A.SiteBasis(), A.Basis1(), B2);

   // loop over components in A and add A \oplus B_zero.
   // then loop over components in B and add A_zero \oplus B
   for (MPOpComponent::const_iterator I = A.begin(); I != A.end(); ++I)
   {
      // make a zero operator over the B basis that transforms as the current component
      MPOpComponent::OperatorType BZero(B.Basis1(), B.Basis2(), I->first);
      Result[I->first] = transform(I->second, ApplyTensorRowSumBind2nd(B2, BZero));
   }

   // do the same for B
   for (MPOpComponent::const_iterator I = B.begin(); I != B.end(); ++I)
   {
      // make a zero operator over the A basis that transforms as the current component
      MPOpComponent::OperatorType AZero(A.Basis1(), A.Basis2(), I->first);
      Result[I->first] += transform(I->second, ApplyTensorRowSumBind1st(B2, AZero));
   }

   Result.debug_check_structure();
   return Result;
}

struct ApplyTensorColSumBind2nd
{
   typedef SimpleOperator const& argument_type;
   typedef SimpleOperator result_type;

   ApplyTensorColSumBind2nd(SumBasis<BasisList> const& B1,
                            SimpleOperator const& Op) : B1_(B1), Op_(Op) {}

   result_type operator()(argument_type x) const
   {
      return tensor_col_sum(x, Op_, B1_);
   }

   SumBasis<BasisList> B1_;
   SimpleOperator Op_;
};

struct ApplyTensorColSumBind1st
{
   typedef SimpleOperator const& argument_type;
   typedef SimpleOperator result_type;

   ApplyTensorColSumBind1st(SumBasis<BasisList> const& B1,
                            SimpleOperator const& Op) : B1_(B1), Op_(Op) {}

   result_type operator()(argument_type x) const
   {
      return tensor_col_sum(Op_, x, B1_);
   }

   SumBasis<BasisList> B1_;
   SimpleOperator Op_;
};

MPOpComponent tensor_col_sum(MPOpComponent const& A,
                             MPOpComponent const& B,
                             SumBasis<BasisList> const& B1)
{
   DEBUG_PRECONDITION_EQUAL(A.SiteBasis(), B.SiteBasis());
   DEBUG_PRECONDITION_EQUAL(A.Basis2(), B.Basis2());

   typedef MPOpComponent::mapped_type OperatorType;

   MPOpComponent Result(A.SiteBasis(), B1, A.Basis2());

   // loop over components in A and add A \oplus B_zero.
   // then loop over components in B and add A_zero \oplus B
   for (MPOpComponent::const_iterator I = A.begin(); I != A.end(); ++I)
   {
      // make a zero operator over the B basis that transforms as the current component
      MPOpComponent::OperatorType BZero(B.Basis1(), B.Basis2(), I->first);
      Result[I->first] = transform(I->second, ApplyTensorColSumBind2nd(B1, BZero));
   }

   // do the same for B
   for (MPOpComponent::const_iterator I = B.begin(); I != B.end(); ++I)
   {
      // make a zero operator over the A basis that transforms as the current component
      MPOpComponent::OperatorType AZero(A.Basis1(), A.Basis2(), I->first);
      Result[I->first] += transform(I->second, ApplyTensorColSumBind1st(B1, AZero));
   }

   Result.debug_check_structure();
   return Result;
}

MPMatrix operator_prod(MPOpComponent const& M,
                       MPStateComponent const& A,
                       MPMatrix const& E,
                       HermitianProxy<MPStateComponent> const& B)
{
   DEBUG_PRECONDITION_EQUAL(M.SiteBasis(), A.SiteBasis());
   DEBUG_PRECONDITION_EQUAL(M.SiteBasis(), B.base().SiteBasis());
   DEBUG_PRECONDITION_EQUAL(M.Basis2(), E.SiteBasis());
   DEBUG_PRECONDITION_EQUAL(A.Basis2(), E.Basis1());
   DEBUG_PRECONDITION_EQUAL(E.Basis2(), B.base().Basis2());

   MPMatrix Result(M.Basis1(), A.Basis1(), B.base().Basis1());
   for (MPOpComponent::const_iterator I = M.begin(); I != M.end(); ++I)
   {
      for (CompoundOperator::const_iterator Ai = iterate(I->second); Ai; ++Ai)
      {
         for (CompoundOperator::const_inner_iterator Aj = iterate(Ai); Aj; ++Aj)
         {

            for (SimpleOperator::const_iterator Aap = iterate(*Aj); Aap; ++Aap)
            {
               for (SimpleOperator::const_inner_iterator Aa = iterate(Aap); Aa; ++Aa)
               {
                  Result[Aa.index1()] += (*Aa) * triple_prod(A[Aj.index1()],
                                                             E[Aa.index2()],
                                                             herm(B.base()[Aj.index2()]),
                                                             I->first,
                                                             M.Basis1()[Aa.index1()]);
               }
            }
         }
      }
   }
   return Result;
}

MPMatrix operator_prod(HermitianProxy<MPOpComponent> const& M,
                       HermitianProxy<MPStateComponent> const& A,
                       MPMatrix const& E,
                       MPStateComponent const& B)
{
   DEBUG_PRECONDITION_EQUAL(M.base().SiteBasis(), A.base().SiteBasis());
   DEBUG_PRECONDITION_EQUAL(M.base().SiteBasis(), B.SiteBasis());
   DEBUG_PRECONDITION_EQUAL(M.base().Basis1(), E.SiteBasis());
   DEBUG_PRECONDITION_EQUAL(A.base().Basis1(), E.Basis1());
   DEBUG_PRECONDITION_EQUAL(E.Basis2(), B.Basis1());

   MPMatrix Result(M.base().Basis2(), A.base().Basis2(), B.Basis2());
   for (MPOpComponent::const_iterator I = M.base().begin(); I != M.base().end(); ++I)
   {
      for (CompoundOperator::const_iterator Ai = iterate(I->second); Ai; ++Ai)
      {
         for (CompoundOperator::const_inner_iterator Aj = iterate(Ai); Aj; ++Aj)
         {

            for (SimpleOperator::const_iterator Aap = iterate(*Aj); Aap; ++Aap)
            {
               for (SimpleOperator::const_inner_iterator Aa = iterate(Aap); Aa; ++Aa)
               {
                  Result[Aa.index2()] += herm(*Aa)
                     * triple_prod(herm(A.base()[Aj.index1()]),
                                   E[Aa.index1()],
                                   B[Aj.index2()],
                                   I->first,
                                   M.base().Basis2()[Aa.index2()]);
               }
            }
         }
      }
   }
   return Result;
}

MPMatrix local_operator_prod(MPOpComponent const& M,
                             MPStateComponent const& E,
                             MPMatrix const& A,
                             HermitianProxy<MPStateComponent> const& F)
{
   DEBUG_PRECONDITION_EQUAL(M.SiteBasis(), A.SiteBasis());
   DEBUG_PRECONDITION_EQUAL(M.Basis1(), E.SiteBasis());
   DEBUG_PRECONDITION_EQUAL(M.Basis2(), F.base().SiteBasis());
   DEBUG_PRECONDITION_EQUAL(E.Basis2(), A.Basis1());
   DEBUG_PRECONDITION_EQUAL(A.Basis2(), F.base().Basis2());

   MPMatrix Result(M.LocalBasis(), E.Basis1(), F.base().Basis1());
   for (MPOpComponent::const_iterator I = M.begin(); I != M.end(); ++I)
   {
      for (CompoundOperator::const_iterator Ai = iterate(I->second); Ai; ++Ai)
      {
         for (CompoundOperator::const_inner_iterator Aj = iterate(Ai); Aj; ++Aj)
         {

            for (SimpleOperator::const_iterator Aap = iterate(*Aj); Aap; ++Aap)
            {
               for (SimpleOperator::const_inner_iterator Aa = iterate(Aap); Aa; ++Aa)
               {
                  Result[Aj.index1()] += (*Aa) * triple_prod(E[Aa.index1()],
                                                             A[Aj.index2()],
                                                             herm(F.base()[Aa.index2()]),
                                                             I->first,
                                                             M.LocalBasis()[Aj.index1()]);
               }
            }
         }
      }
   }
   return Result;
}

#if 0
MPMatrix local_operator_prod(HermitianProxy<MPOpComponent> const& M,
                             HermitianProxy<MPStateComponent> const& E,
                             MPMatrix const& A,
                             MPStateComponent const& F)
{
   DEBUG_PRECONDITION_EQUAL(M.base().SiteBasis(), A.base().SiteBasis());
   DEBUG_PRECONDITION_EQUAL(M.base().SiteBasis(), B.SiteBasis());
   DEBUG_PRECONDITION_EQUAL(M.base().Basis1(), E.SiteBasis());
   DEBUG_PRECONDITION_EQUAL(A.base().Basis1(), E.Basis1());
   DEBUG_PRECONDITION_EQUAL(E.Basis2(), B.Basis1());

   MPMatrix Result(M.base().Basis2(), A.base().Basis2(), B.Basis2());
   for (MPOpComponent::const_iterator I = M.base().begin(); I != M.base().end(); ++I)
   {
      for (CompoundOperator::const_iterator Ai = iterate(I->second); Ai; ++Ai)
      {
         for (CompoundOperator::const_inner_iterator Aj = iterate(Ai); Aj; ++Aj)
         {

            for (SimpleOperator::const_iterator Aap = iterate(*Aj); Aap; ++Aap)
            {
               for (SimpleOperator::const_inner_iterator Aa = iterate(Aap); Aa; ++Aa)
               {
                  Result[Aa.index2()] += herm(*Aa)
                     * triple_prod(herm(A.base()[Aj.index1()]),
                                   E[Aa.index1()],
                                   B[Aj.index2()],
                                   I->first,
                                   M.base().Basis2()[Aa.index2()]);
               }
            }
         }
      }
   }
   return Result;
}
#endif

MPOpComponent prod(MPOpComponent const& A, SimpleOperator const& Op)
{
   PRECONDITION(is_scalar(Op.TransformsAs()))(Op.TransformsAs());
   DEBUG_PRECONDITION_EQUAL(A.Basis2(), Op.Basis1());

   MPOpComponent Result(A.SiteBasis(), A.Basis1(), Op.Basis2());
   for (MPOpComponent::const_iterator AI = A.begin(); AI != A.end(); ++AI)
   {
      Result[AI->first] = transform(AI->second,
         bind_second(LinearAlgebra::IrredProd<SimpleOperator, SimpleOperator>(AI->first), Op));

      //prod(AI->second, Op, AI->first);
   }

   Result.debug_check_structure();
   return Result;
}

MPOpComponent prod(SimpleOperator const& Op, MPOpComponent const& A)
{
   A.debug_check_structure();

   PRECONDITION(is_scalar(Op.TransformsAs()))(Op.TransformsAs());
   DEBUG_PRECONDITION_EQUAL(Op.Basis2(), A.Basis1());

   MPOpComponent Result(A.SiteBasis(), Op.Basis1(), A.Basis2());
   for (MPOpComponent::const_iterator AI = A.begin(); AI != A.end(); ++AI)
   {
      Result[AI->first] = transform(AI->second,
         bind_first(LinearAlgebra::IrredProd<SimpleOperator, SimpleOperator>(AI->first), Op));
   }

   Result.debug_check_structure();
   return Result;
}

MPOpComponent prod(MPOpComponent const& A, HermitianProxy<SimpleOperator> const& Op)
{
   PRECONDITION(is_scalar(Op.base().TransformsAs()))(Op.base().TransformsAs());
   DEBUG_PRECONDITION_EQUAL(A.Basis2(), Op.base().Basis2());

   MPOpComponent Result(A.SiteBasis(), A.Basis1(), Op.base().Basis1());
   for (MPOpComponent::const_iterator AI = A.begin(); AI != A.end(); ++AI)
   {
      Result[AI->first] = transform(AI->second,
         bind_second(IrredProd_Herm<SimpleOperator, HermitianProxy<SimpleOperator> >(), Op));

      //prod(AI->second, Op, AI->first);
   }

   Result.debug_check_structure();
   return Result;
}

MPOpComponent prod(HermitianProxy<SimpleOperator> const& Op, MPOpComponent const& A)
{
   A.debug_check_structure();

   PRECONDITION(is_scalar(Op.base().TransformsAs()))(Op.base().TransformsAs());
   DEBUG_PRECONDITION_EQUAL(Op.base().Basis1(), A.Basis1());

   MPOpComponent Result(A.SiteBasis(), Op.base().Basis2(), A.Basis2());
   for (MPOpComponent::const_iterator AI = A.begin(); AI != A.end(); ++AI)
   {
      Result[AI->first] = transform(AI->second,
         bind_first(IrredProd_Herm<HermitianProxy<SimpleOperator>, SimpleOperator>(), Op));
   }

   Result.debug_check_structure();
   return Result;
}

MPOpComponent local_tensor_prod(MPOpComponent const& A, MPOpComponent const& B)
{
   using QuantumNumbers::QuantumNumberList;
   typedef LinearAlgebra::IrredProd<SimpleOperator, SimpleOperator> IrredProdSimpleOperator;
   Tensor::ProductBasis<BasisList, BasisList> PB(A.LocalBasis(), B.LocalBasis());
   MPOpComponent Result(PB.Basis(), A.Basis1(), B.Basis2());
   for (MPOpComponent::const_iterator AI = A.begin(); AI != A.end(); ++AI)
   {
      for (MPOpComponent::const_iterator BI = B.begin(); BI != B.end(); ++BI)
      {
         QuantumNumberList ql = transform_targets(AI->first, BI->first);
         for (QuantumNumberList::const_iterator qI = ql.begin(); qI != ql.end(); ++qI)
         {
            Result[*qI] += tensor_prod(AI->second, BI->second, PB, PB, *qI,
                                       IrredProdSimpleOperator(*qI));
         }
      }
   }
   return Result;
}

SimpleOperator IdentityOperator(BasisList const& B)
{
   SimpleOperator Result(B, B, QuantumNumbers::QuantumNumber(B.GetSymmetryList()));
   for (std::size_t i = 0; i < B.size(); ++i)
   {
      Result(i,i) = 1.0;
   }
   return Result;
}

bool IsProportionalIdentity(MPOpComponent const& Op)
{
   if (Op.Basis1().size() != 1 || !is_scalar(Op.Basis1()[0])
       || Op.Basis2().size() != 1 || !is_scalar(Op.Basis2()[0])) return false;

   MPOpComponent Ident = ConstructIdentity(Op.SiteBasis());

   double Scale = trace(scalar_prod(Ident, herm(Ident))).real();
   double Factor = trace(scalar_prod(Op, herm(Op))).real();
   std::complex<double> Fac = trace(scalar_prod(Op, herm(Ident)));
   double Fac2 = std::norm(Fac) / Scale;

   //TRACE(Factor)(Fac2)(Scale)(Fac)(Fac/Scale);

   return (std::abs(Factor - Fac2)
           < 100 * std::numeric_limits<double>::epsilon() * (1 + std::abs(Factor+Fac2)));
}

std::complex<double>
IdentityScale(MPOpComponent const& Op)
{
   DEBUG_CHECK(IsProportionalIdentity(Op));

   MPOpComponent Ident = ConstructIdentity(Op.SiteBasis());
   return trace(scalar_prod(Op, herm(Ident))) / trace(scalar_prod(Ident, herm(Ident))).real();
}

SimpleOperator ExpandBasis1(MPOpComponent& A, Normalization n)
{
   //   PANIC("*******WARNING*********** ExpandBasis1");
   // iterate through A and determine which basis states need to be in the full basis
   BasisList LocalBasis(A.GetSymmetryList());
   std::vector<std::pair<int, int> > LocalRMap;
   for (MPOpComponent::const_iterator MI = A.begin(); MI != A.end(); ++MI)
   {
      for (LinearAlgebra::const_iterator<CompoundOperator>::type I = iterate(MI->second); I; ++I)
      {
         for (LinearAlgebra::const_inner_iterator<CompoundOperator>::type J = iterate(I); J; ++J)
         {
            LocalBasis.push_back(MI->first);
            LocalRMap.push_back(std::make_pair(J.index1(), J.index2()));
         }
      }
   }

   ProductBasis<BasisList, BasisList> FullBasis1(LocalBasis, A.Basis2());
   MPOpComponent Result(A.SiteBasis(), FullBasis1.Basis(), A.Basis2());

   for (std::size_t t = 0; t < FullBasis1.size(); ++t)
   {
      int l, b2;
      std::tie(l,b2) = FullBasis1.rmap(t);

      int sp, s;
      std::tie(sp, s) = LocalRMap[l];

      if (Result[LocalBasis[l]](sp, s).is_null())
      {
         Result[LocalBasis[l]](sp, s) = SimpleOperator(FullBasis1.Basis(), A.Basis2(), LocalBasis[l]);
      }

      Result[LocalBasis[l]](sp, s)(t, b2) = 1.0 / std::sqrt(double(degree(A.SiteBasis()[sp])));
   }

   // check the normalization
   DEBUG_CHECK_CLOSE(norm_frob_sq(scalar_prod(Result, herm(Result))),
                     VectorBasis(FullBasis1.Basis()).total_degree());

   //   TRACE(A);
   //   TRACE(scalar_prod(Result, herm(Result)));

   SimpleOperator Res = scalar_prod(A, herm(Result));

   //   double RDeg = std::sqrt(double(A.SiteBasis().total_degree()));
   //   Result *= RDeg;
   //   Res *= 1.0 / RDeg;

   //   Result = prod(Res, Result);
   //   Res = IdentityOperator(Res.Basis1());

   //   TRACE(Result)(Res);

   //   TRACE(Res)(prod(Res, Result));
   A = Result;
   return Res;
}

SimpleOperator ExpandBasis2(MPOpComponent& A, Normalization n)
{
   //   PANIC("*******WARNING*********** ExpandBasis2");
   // iterate through A and determine which basis states need to be in the full basis
   BasisList LocalBasis(A.GetSymmetryList());
   std::vector<std::pair<int, int> > LocalRMap;
   for (MPOpComponent::const_iterator MI = A.begin(); MI != A.end(); ++MI)
   {
      for (LinearAlgebra::const_iterator<CompoundOperator>::type I = iterate(MI->second); I; ++I)
      {
         for (LinearAlgebra::const_inner_iterator<CompoundOperator>::type J = iterate(I); J; ++J)
         {
            LocalBasis.push_back(MI->first);
            LocalRMap.push_back(std::make_pair(J.index1(), J.index2()));
         }
      }
   }

   ProductBasis<BasisList, BasisList> FullBasis2(A.Basis1(), adjoint(LocalBasis));
   MPOpComponent Result(A.SiteBasis(), A.Basis1(), FullBasis2.Basis());

   for (std::size_t t = 0; t < FullBasis2.size(); ++t)
   {
      int l, b1;
      std::tie(b1,l) = FullBasis2.rmap(t);

      int sp, s;
      std::tie(sp, s) = LocalRMap[l];

      if (Result[LocalBasis[l]](sp, s).is_null())
      {
         Result[LocalBasis[l]](sp, s)
            = SimpleOperator(A.Basis1(), FullBasis2.Basis(), LocalBasis[l]);
      }

      Result[LocalBasis[l]](sp, s)(b1, t) =
         std::sqrt(double(degree(FullBasis2[t]))
                   / (degree(A.Basis1()[b1]) * degree(A.SiteBasis()[sp])));
   }

   //   TRACE(A);
   //   TRACE(scalar_prod(herm(Result), Result));

   // check the normalization
   DEBUG_CHECK_CLOSE(norm_frob_sq(scalar_prod(herm(Result), Result)),
                     VectorBasis(FullBasis2.Basis()).total_degree());

   SimpleOperator Res = scalar_prod(herm(Result), A);

   //   Result = prod(Result, Res);
   //   Res = IdentityOperator(Res.Basis2());

   //   double RDeg = std::sqrt(double(A.SiteBasis().total_degree()));
   //   Result *= RDeg;
   //   Res *= 1.0 / RDeg;

   //   TRACE(Result)(Res);

   //   TRACE(Res)(prod(Result, Res));
   A = Result;
   return Res;
}

SimpleOperator
local_trace(MPOpComponent const& A)
{
   SimpleOperator Result(A.Basis1(), A.Basis2(), QuantumNumbers::QuantumNumber(A.GetSymmetryList()));
   MPOpComponent::mapped_type x = A[QuantumNumbers::QuantumNumber(A.GetSymmetryList())];
   for (MPOpComponent::mapped_type::const_iterator I = iterate(x); I; ++I)
   {
      for (MPOpComponent::mapped_type::const_inner_iterator J = iterate(I); J; ++J)
      {
         if (J.index1() == J.index2())
            Result += trace(x.Basis1()[J.index1()]) * (*J);
      }
   }
   return Result;
}

double const OverlapEpsilon = 1E-14;       // relative epsilon for parallel entries
double const SmallEpsilon = 1E-14;  // absolute epsilon for zero entries

MPOpComponent
SumFix1(MPOpComponent const& A, MPOpComponent const& B,
        SimpleOperator& AMap, SimpleOperator& BMap)
{
   DEBUG_CHECK_EQUAL(A.Basis1(), B.Basis1());

   BasisList NewBasis2(A.Basis2());
   for (unsigned i = 0; i < B.Basis2().size(); ++i)
   {
      NewBasis2.push_back(B.Basis2()[i]);
   }

   AMap = SimpleOperator(NewBasis2, A.Basis2());
   for (unsigned i = 0; i < A.Basis2().size(); ++i)
   {
      AMap(i,i) = 1.0;
   }

   int ASize = A.Basis2().size();
   BMap = SimpleOperator(NewBasis2, B.Basis2());
   for (unsigned i = 0; i < B.Basis2().size(); ++i)
   {
      BMap(i+ASize,i) = 1.0;
   }

   MPOpComponent Result = A*herm(AMap) + B*herm(BMap);
   SimpleOperator Res = TruncateBasis2(Result);
   AMap = Res*AMap;
   BMap = Res*BMap;
   return Result;
}

SimpleOperator TruncateBasis1(MPOpComponent& A)
{
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
      A = MPOpComponent(A.SiteBasis(), A.Basis2(), A.Basis2());
      return Result;
   }

   // Normalize the matrices by the degree of the site basis
   // This keeps all matrix elements ~ 1
   double const ASiteBasisDegree = A.SiteBasis().total_degree();

   SimpleOperator Overlaps = scalar_prod(A, herm(A));

   // Normalize each row
   SimpleOperator Normalizer(A.Basis1(), A.Basis1());
   SimpleOperator InverseNormalizer(A.Basis1(), A.Basis1());
   for (unsigned i = 0; i < A.Basis1().size(); ++i)
   {
      double x = std::sqrt(Overlaps(i,i).real() / ASiteBasisDegree);
      if (x <= OperatorTruncationEpsilon)
      {
         // one of the components is zero.  No need to remove it here,
         // it will get removed later anyway.
         x = 1;
      }
      InverseNormalizer(i,i) = x;
      Normalizer(i,i) = 1.0 / x;
   }
   MPOpComponent AOld = A;
   A = prod(Normalizer, A);
   Overlaps = scalar_prod(A, herm(A));

   // This is the transform that truncates the rows of A.
   // row[i] = NewRows[i].second * A(NewRows[i].second, all)
   std::vector<std::pair<int, std::complex<double> > > NewRows;
   std::set<int> KeepStates; // the states that we are going to keep
   for (int i = 0; i < int(A.Basis1().size()); ++i)
   {
      NewRows.push_back(std::make_pair(i, 1.0));
      double imat = Overlaps(i,i).real();
      // if the row is zero, we can eliminate it completely.
      // Because we've normalized everything, then it is either 1 or ~epsilon here.
      if (imat <= 0.1)
      {
         NewRows.back().first = -1;  // special value, indicates the row is not needed
         continue;
      }
      bool Parallel = false;  // is vector i parallel to some other vector?
      // loop to find out if row i is parallel to row j
      for (int j = 0; j < i; ++j)
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

   // adjust for the normalizer
   Reg = InverseNormalizer * Reg;

   MPOpComponent tA = prod(Trunc, A); // the result

   // remove any components tA[n] that are empty
   MPOpComponent::iterator I = tA.begin();
   while (I != tA.end())
   {
      if (norm_frob_sq(I->second) == 0)
         tA.Data_.erase(I++);
      else
         ++I;
   }

#if !defined(NDEBUG)
   // verify that prod(Reg, tA) is the same as A.  This is slightly awkward
   // as we don't have operator-() for MPOpComponent
   MPOpComponent ACheck = prod(Reg, tA);
   std::complex<double> a = trace(scalar_prod(ACheck, herm(ACheck)));
   std::complex<double> b = trace(scalar_prod(AOld, herm(AOld)));
   std::complex<double> c = trace(scalar_prod(AOld, herm(ACheck)));
   CHECK(LinearAlgebra::norm_frob(a.real()+b.real()-2.0*c.real()) <= 1E-10 * (norm_frob(a.real()+norm_frob(b.real()))))
      (a)(b)(c)(A)(AOld)
      (tA)(ACheck)(Trunc)(Reg)(Overlaps);
#endif

   A = tA;
   return Reg;
}

SimpleOperator TruncateBasis2(MPOpComponent& A)
{
   // if the operator is trivially zero, then return early
   if (A.Basis2().size() == 0)
   {
      return SimpleOperator(A.Basis2(), A.Basis2(), QuantumNumber(A.GetSymmetryList()));
   }
   else if (A.Basis1().size() == 0)
   {
      SimpleOperator Result(A.Basis1(), A.Basis2(), QuantumNumber(A.GetSymmetryList()));
      A = MPOpComponent(A.SiteBasis(), A.Basis1(), A.Basis1());
      return Result;
   }

   // Normalize the matrices by the degree of the site basis.
   // This keeps all matrix elements ~ 1
   double const ASiteBasisDegree = A.SiteBasis().total_degree();

   SimpleOperator Overlaps = scalar_prod(herm(A), A);

   // Normalize each column
   SimpleOperator Normalizer(A.Basis2(), A.Basis2());
   SimpleOperator InverseNormalizer(A.Basis2(), A.Basis2());
   for (unsigned i = 0; i < A.Basis2().size(); ++i)
   {
      double x = std::sqrt(Overlaps(i,i).real() / ASiteBasisDegree);
      if (x <= OperatorTruncationEpsilon)
      {
         // one of the components is zero.  No need to remove it here,
         // it will get removed later anyway.
         x = 1;
      }
      InverseNormalizer(i,i) = x;
      Normalizer(i,i) = 1.0 / x;
   }
   MPOpComponent AOld = A;
   A = prod(A, Normalizer);
   Overlaps = scalar_prod(herm(A), A);

   // This is the transform that truncates the columns of A.
   // row[i] = NewCols[i].second * A(all, NewCols[i].second)
   std::vector<std::pair<int, std::complex<double> > > NewCols;
   std::set<int> KeepStates; // the states that we are going to keep
   for (int i = 0; i < int(A.Basis2().size()); ++i)
   {
      NewCols.push_back(std::make_pair(i, 1.0));
      double imat = Overlaps(i,i).real();
      // if the row is zero, we can eliminate it completely
      // Because we've normalized everything, then it is either 1 or ~epsilon here.
      if (imat <= 0.1)
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
   Reg = Reg * InverseNormalizer;

   MPOpComponent tA = prod(A, Trunc); // the result

   // remove any components tA[n] that are empty
   MPOpComponent::iterator I = tA.begin();
   while (I != tA.end())
   {
      if (norm_frob_sq(I->second) == 0)
         tA.Data_.erase(I++);
      else
         ++I;
   }

#if !defined(NDEBUG)
   // verify that prod(tA, Reg) is the same as A.  This is slightly awkward
   // as we don't have operator-() for MPOpComponent
   MPOpComponent ACheck = prod(tA, Reg);
   std::complex<double> a = trace(scalar_prod(herm(ACheck), ACheck));
   std::complex<double> b = trace(scalar_prod(herm(AOld), AOld));
   std::complex<double> c = trace(scalar_prod(herm(AOld), ACheck));
   CHECK(LinearAlgebra::norm_frob(a.real()+b.real()-2.0*c.real()) <= 1E-10 * (norm_frob(a.real()+norm_frob(b.real()))))
      (a)(b)(c)(A)(AOld)
      (tA)(ACheck)(Trunc)(Reg)(Overlaps);
#endif

   A = tA;
   return Reg;
}


MPStateComponent
MProd<MPOpComponent, MPStateComponent>::operator()(MPOpComponent const& M,
                                                   MPStateComponent const& A,
                                                   ProductBasis<BasisList, VectorBasis> const& B1,
                                                   ProductBasis<BasisList, VectorBasis> const& B2)
{
   MPStateComponent Result(M.SiteBasis(), B1.Basis(), B2.Basis());

   for (MPOpComponent::const_iterator MI = M.begin(); MI != M.end(); ++MI)
   {
      for (const_iterator<CompoundOperator>::type I = iterate(MI->second); I; ++I)
      {
         for (const_inner_iterator<CompoundOperator>::type J = iterate(I); J; ++J)
         {
            Result[J.index1()] += QuantumNumbers::conj_phase(M.SiteBasis()[J.index2()],
                                                             MI->first,
                                                             M.SiteBasis()[J.index1()])
               * Tensor::tensor_prod(*J, A[J.index2()], B1, B2, M.SiteBasis()[J.index1()]);
         }
      }
   }
   return Result;
}

MPOpComponent
MProd<MPOpComponent, MPOpComponent>::operator()(MPOpComponent const& M,
                                                MPOpComponent const& N,
                                                ProductBasis<BasisList, BasisList> const& B1,
                                                ProductBasis<BasisList, BasisList> const& B2)
{
   MPOpComponent Result(M.SiteBasis(), B1.Basis(), B2.Basis());

   for (MPOpComponent::const_iterator MI = M.begin(); MI != M.end(); ++MI)
   {
      for (MPOpComponent::const_iterator NI = N.begin(); NI != N.end(); ++NI)
      {
         QuantumNumbers::QuantumNumberList QL(transform_targets(MI->first, NI->first));
         for (std::size_t i = 0; i < QL.size(); ++i)
         {
            Result[QL[i]] += prod(MI->second, NI->second, QL[i],
                                  Tensor::TensorProd<SimpleOperator, SimpleOperator>(B1, B2, QL[i]));
         }
      }
   }
   return Result;
}

MPOpComponent conj(MPOpComponent const& x)
{
   MPOpComponent Result(x);
   for (MPOpComponent::iterator I = Result.begin(); I != Result.end(); ++I)
   {
      I->second = conj(I->second);
   }
   return Result;
}

MPOpComponent local_adjoint(MPOpComponent const& x)
{
   MPOpComponent Result(x.SiteBasis(), adjoint(x.Basis1()), adjoint(x.Basis2()));
   for (MPOpComponent::const_iterator I = x.begin(); I != x.end(); ++I)
   {
      Result[adjoint(I->first)] = adjoint(I->second, LinearAlgebra::FlipConjugate<MPOpComponent::OperatorType>());
   }
   return Result;
}

LinearAlgebra::Vector<SimpleOperator>
project_row(MPOpComponent const& x, int r)
{
   LinearAlgebra::Vector<SimpleOperator> Result(x.Basis2().size());
   for (MPOpComponent::const_iterator MI = x.begin(); MI != x.end(); ++MI)
   {
      for (const_iterator<CompoundOperator>::type I = iterate(MI->second); I; ++I)
      {
         for (const_inner_iterator<CompoundOperator>::type J = iterate(I); J; ++J)
         {
            for (const_iterator<SimpleOperator>::type K = iterate(*J); K; ++K)
            {
               for (const_inner_iterator<SimpleOperator>::type L = iterate(K); L; ++L)
               {
                  if (int(L.index1()) == r)
                  {
                     if (Result[L.index2()].is_null())
                     {
                        Result[L.index2()] = SimpleOperator(x.LocalBasis(), x.LocalBasis(), MI->first);
                     }
                     else
                     {
                        CHECK_EQUAL(Result[L.index2()].TransformsAs(), MI->first)
                           ("MPOpComponent is not irreducible");
                     }
                     Result[L.index2()](J.index1(), J.index2()) = *L;
                  }
               }
            }
         }
      }
   }
   return Result;
}

LinearAlgebra::Vector<SimpleOperator>
project_column(MPOpComponent const& x, int c)
{
   LinearAlgebra::Vector<SimpleOperator> Result(x.Basis1().size());
   for (MPOpComponent::const_iterator MI = x.begin(); MI != x.end(); ++MI)
   {
      for (const_iterator<CompoundOperator>::type I = iterate(MI->second); I; ++I)
      {
         for (const_inner_iterator<CompoundOperator>::type J = iterate(I); J; ++J)
         {
            for (const_iterator<SimpleOperator>::type K = iterate(*J); K; ++K)
            {
               for (const_inner_iterator<SimpleOperator>::type L = iterate(K); L; ++L)
               {
                  if (int(L.index2()) == c)
                  {
                     if (Result[L.index1()].is_null())
                     {
                        Result[L.index1()] = SimpleOperator(x.LocalBasis(), x.LocalBasis(), MI->first);
                     }
                     else
                     {
                        CHECK_EQUAL(Result[L.index1()].TransformsAs(), MI->first)
                           ("MPOpComponent is not irreducible");
                     }
                     Result[L.index1()](J.index1(), J.index2()) = *L;
                  }
               }
            }
         }
      }
   }
   return Result;
}

#if 0
MPOpMatrix as_matrix(MPOpComponent::ComponentType const& x)
{
   // get the first non-zero element so we know what the matrix basis is
   const_iterator<MPOpComponent::ComponentType>::type I = iterate(x);
   if (!I) return MPOpMatrix();
   const_inner_iterator<MPOpComponent::ComponentType>::type J = iterate(I);
   if (!J) return MPOpMatrix();
   MPOpMatrix Result(J->Basis1(), J->Basis2(), x.TransformsAs());

   for ( ; I; ++I)
   {
      for (J = iterate(I); J; ++J)
      {
         for (const_iterator<SimpleOperator>::type K = iterate(*J); K; ++K)
            {
               for (const_inner_iterator<SimpleOperator>::type L = iterate(K); L; ++L)
               {
                  if (Result(L.index1(), L.index2()).is_null())
                     Result(L.index1(), L.index2()) = SimpleOperator(x.Basis1(), x.Basis2(), x.TransformsAs());
                  Result(L.index1(), L.index2())(J.index1(), J.index2()) = *L;
               }
            }
      }
   }
   return Result;
}
#endif
