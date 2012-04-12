// -*- C++ -*- $Id$

#include "mpopcomponent.h"
#include <boost/tuple/tuple.hpp>
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
