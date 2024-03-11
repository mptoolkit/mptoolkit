// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mpo/generic_mpo.cpp
//
// Copyright (C) 2012-2022 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2022 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#include "generic_mpo.h"

std::vector<BasisList>
GenericMPO::LocalBasis1List() const
{
   std::vector<BasisList> Result;
   Result.reserve(Data_.size());
   for (unsigned i = 0; i < Data_.size(); ++i)
   {
      Result.push_back(Data_[i].LocalBasis1());
   }
   return Result;
}

std::vector<BasisList>
GenericMPO::LocalBasis2List() const
{
   std::vector<BasisList> Result;
   Result.reserve(Data_.size());
   for (unsigned i = 0; i < Data_.size(); ++i)
   {
      Result.push_back(Data_[i].LocalBasis2());
   }
   return Result;
}

void
GenericMPO::check_structure() const
{
   const_iterator I = this->begin();
   while (I != this->end())
   {
      I->check_structure();
      const_iterator J = I;
      ++I;
      if (I != this->end())
      {
         CHECK_EQUAL(J->Basis2(), I->Basis1());
      }
   }
}

std::ostream&
operator<<(std::ostream& out, GenericMPO const& op)
{
   out << "Operator has a unit cell of " << op.size() << " sites.\n";
   for (unsigned i = 0; i < op.size(); ++i)
   {
      out << "Site " << i << ": " << op[i] << '\n';
   }
   return out;
}

bool
GenericMPO::is_null() const
{
   for (unsigned i = 0; i < Data_.size(); ++i)
   {
      if (!Data_[i].is_null())
         return false;
   }
   return true;
}

PStream::opstream& operator<<(PStream::opstream& out, GenericMPO const& op)
{
   return out << op.Data_;
}

PStream::ipstream& operator>>(PStream::ipstream& in, GenericMPO& op)
{
   return in >> op.Data_;
}

#if 0
GenericMPO&
operator*=(GenericMPO& x, double a)
{
   x.front() *= a;
   return x;
}

GenericMPO&
operator*=(GenericMPO& x, std::complex<double> a)
{
   x.front() *= a;
   return x;
}

GenericMPO operator*(double a, GenericMPO const& x)
{
   GenericMPO Result(x);
   Result *= a;
   return Result;
}

GenericMPO operator*(GenericMPO const& x, double a)
{
   GenericMPO Result(x);
   Result *= a;
   return Result;
}

GenericMPO operator*(std::complex<double> a, GenericMPO const& x)
{
   GenericMPO Result(x);
   Result *= a;
   return Result;
}

GenericMPO operator*(GenericMPO const& x, std::complex<double> a)
{
   GenericMPO Result(x);
   Result *= a;
   return Result;
}
#endif

void zero_unused_elements(GenericMPO& Op)
{
   bool Done = false;
   while (!Done)
   {
      Done = true;
      std::set<int> NextKeep;
      GenericMPO::iterator I = Op.begin();

      std::set<int> RowsToKeep;
      for (unsigned i = 0; i < I->size1(); ++i)
         for (unsigned j = 0; j < I->size2(); ++j)
            if (I->iterate_at(i,j))
               {
                  RowsToKeep.insert(j);
               }
      ++I;

      while (I != Op.end())
      {
         OperatorComponent Current = *I;
         for (unsigned i = 0; i < Current.size1(); ++i)
         {
            for (unsigned j = 0; j < Current.size2(); ++j)
            {
               if (!Current.iterate_at(i,j))
                  continue;

               if (RowsToKeep.count(i))
                  NextKeep.insert(j);
               else
               {
                  zero_element(Current.data(), i,j);
                  Done = false;
               }

            }
         }
         *I = Current;

         RowsToKeep = NextKeep;
         NextKeep.clear();

         ++I;
      }

      // now work backwards
      --I;
      std::set<int> ColumnsToKeep;
      for (unsigned i = 0; i < I->size1(); ++i)
         for (unsigned j = 0; j < I->size2(); ++j)
            if (I->iterate_at(i,j))
               ColumnsToKeep.insert(i);

      while (I != Op.begin())
      {
         --I;

         OperatorComponent Current = *I;
         for (unsigned i = 0; i < Current.size1(); ++i)
         {
            for (unsigned j = 0; j < Current.size2(); ++j)
            {
               if (!Current.iterate_at(i,j))
                  continue;

               if (ColumnsToKeep.count(j))
                  NextKeep.insert(i);
               else
               {
                  zero_element(Current.data(), i,j);
                  Done = false;
               }

            }
         }
         *I = Current;

         ColumnsToKeep = NextKeep;
         NextKeep.clear();
      }

   } // while (!Done)
}

bool update_mask(OperatorComponent const& x, OperatorComponent const& y,
                 std::vector<int> const& M1, std::vector<int>& Mask, std::vector<int> const& M2)
{
   typedef OperatorComponent::data_type OpMatrixType;
   typedef OpMatrixType::outer_value_type VecOfVecType;
   //   typedef VecOfVecType::value_type InnerType;
   typedef LinearAlgebra::MapVector<SimpleRedOperator> InnerType;

   bool Updated = false;
   for (const_iterator<OperatorComponent::data_type>::type I = iterate(x.data()); I; ++I)
   {
      if (!M1[I.index()])
         continue;

      for (const_inner_iterator<OperatorComponent::data_type>::type J = iterate(I); J; ++J)
      {
         //         if (M1[J.index1()] && Mask[J.index2()])
         if (Mask[J.index2()])
         {
            const_iterator<InnerType>::type K = iterate(y.data().vec()[J.index2()]);
            while (K && !M2[K.index()])
               ++K;

            if (!K)
            {
               Mask[J.index2()] = 0;
               Updated = true;
            }
         }
      }
   }

   return Updated;
}

bool update_mask(OperatorComponent const& x, std::vector<int>& M1, std::vector<int>& M2)
{
   typedef OperatorComponent::data_type OpMatrixType;
   typedef OpMatrixType::outer_value_type VecOfVecType;
   //   typedef VecOfVecType::value_type InnerType;
   typedef LinearAlgebra::MapVector<SimpleRedOperator> InnerType;

   bool Updated = false;
   for (const_iterator<OperatorComponent::data_type>::type I = iterate(x.data()); I; ++I)
   {
      if (!M1[I.index()])
         continue;

      bool Found = false;
      for (const_inner_iterator<OperatorComponent::data_type>::type J = iterate(I); J; ++J)
      {
         if (M2[J.index2()])
         {
            M2[J.index2()] = 2;  // tag as accessible
            Found = true;
         }
      }
      if (!Found)
      {
         M1[I.index()] = 0;
         Updated = true;
      }
   }

   for (unsigned j = 0; j < M2.size(); ++j)
   {
      if (M2[j] == 2)
      {
         M2[j] = 1;
      }
      else if (M2[j] == 1)
      {
         M2[j] = 0;  // not accessible
         Updated = true;
      }
   }

   return Updated;
}

// remove unused components from the shared basis between x and y
bool cull_boundary(OperatorComponent& x, OperatorComponent& y)
{
   std::set<int> ToKeep;
   for (const_iterator<OperatorComponent::data_type>::type I = iterate(x.data()); I; ++I)
   {
      for (const_inner_iterator<OperatorComponent::data_type>::type J = iterate(I); J; ++J)
      {
         if (nnz(y.data().vec()[J.index2()]) > 0)
            ToKeep.insert(J.index2());
      }
   }
   if (ToKeep.size() != x.Basis2().size())
   {
      x = project_columns(x, ToKeep);
      y = project_rows(y, ToKeep);
      return true;
   }
   return false;
}

void cull_unused_elements(GenericMPO& Op)
{
   // We need at least one bond to optimize
   if (Op.size() < 2)
      return;

   bool Done = false;
   while (!Done)
   {
      Done = true;
      GenericMPO::iterator I = Op.begin();
      GenericMPO::iterator J = I; ++J;

      while (J != Op.end())
      {
         if (cull_boundary(*I, *J))
            Done = false;

         I=J;
         ++J;
      }

      // now work backwards and cull columns
      while (I != Op.begin())
      {
         J=I;
         --I;

         if (cull_boundary(*I, *J))
            Done = false;
      }
   } // while (!Done)
}

void mask_unused_elements(GenericMPO const& Op, std::vector<std::vector<int> >& Mask)
{
   if (Op.size() < 1)
      return;

   GenericMPO::const_iterator I = Op.begin();
   std::vector<std::vector<int> >::iterator M1 = Mask.begin();
   std::vector<std::vector<int> >::iterator M2 = Mask.begin()+1;

   bool Done = false;
   while (!Done)
   {
      Done = true;
      while (I != Op.end())
      {
         if (update_mask(*I, *M1, *M2))
            Done = false;

         ++I;
         ++M1;
         ++M2;
      }

      if (Done)
         break;

      while (I != Op.begin())
      {
         --I;
         --M1;
         --M2;

         if (update_mask(*I, *M1, *M2))
            Done = false;
      }
   }
}

void initialize_mask(GenericMPO const& Op, std::vector<std::vector<int> >& Mask)
{
   std::vector<std::vector<int> >(Op.size()+1).swap(Mask);
   for (unsigned i = 0; i < Op.size(); ++i)
   {
      Mask[i] = std::vector<int>(Op[i].Basis1().size(), true);
   }
   Mask.back() = std::vector<int>(Op.Basis2().size(), true);
}

std::vector<std::vector<int> >
mask_column(GenericMPO const& Op, int Col)
{
   std::vector<std::vector<int> > Mask;
   initialize_mask(Op, Mask);
   std::fill(Mask.back().begin(), Mask.back().end(), false);
   Mask.back()[Col] = true;
   mask_unused_elements(Op, Mask);
   return Mask;
}

std::vector<std::vector<int> >
mask_row(GenericMPO const& Op, int Row)
{
   std::vector<std::vector<int> > Mask;
   initialize_mask(Op, Mask);
   std::fill(Mask.front().begin(), Mask.front().end(), false);
   Mask.front()[Row] = true;
   mask_unused_elements(Op, Mask);
   return Mask;
}

GenericMPO coarse_grain(GenericMPO const& Op, int N)
{
   if (N == 1)
      return Op;

   CHECK(Op.size() % N == 0);
   std::vector<OperatorComponent> Result;
   for (unsigned i = 0; i < Op.size(); i += N)
   {
      OperatorComponent c = local_tensor_prod(Op[i], Op[i+1]);
      for (int j = 2; j < N; ++j)
      {
         c = local_tensor_prod(c, Op[i+j]);
      }
      Result.push_back(c);
   }
   return GenericMPO(Result.begin(), Result.end());
}

GenericMPO coarse_grain_range(GenericMPO const& Op, int beg, int end)
{
   CHECK(0 <= beg && beg <= end && end <= int(Op.size()));
   // coarse-graining 0 or 1 sites is a no-op
   if (end-beg < 2)
      return Op;

   // sites up to beg are unaffected
   std::vector<OperatorComponent> Result;
   for (int i = 0; i < beg; ++i)
   {
      Result.push_back(Op[i]);
   }

   // coarse-grain [beg,end)
   OperatorComponent c = local_tensor_prod(Op[beg], Op[beg+1]);
   for (int j = beg+2; j < end; ++j)
   {
      c = local_tensor_prod(c, Op[j]);
   }
   Result.push_back(c);

   // sites [end, Op.size) are unaffected
   for (int i = end; i < int(Op.size()); ++i)
   {
      Result.push_back(Op[i]);
   }
   return GenericMPO(Result.begin(), Result.end());
}

SimpleOperator make_projector_onto(BasisList const& Basis, std::set<int> const& Onto)
{
   BasisList ProjectedBasis(Basis.GetSymmetryList());
   for (std::set<int>::const_iterator I = Onto.begin(); I != Onto.end(); ++I)
   {
      ProjectedBasis.push_back(Basis[*I]);
   }

   SimpleOperator Result(ProjectedBasis, Basis);
   int j = 0;
   for (std::set<int>::const_iterator I = Onto.begin(); I != Onto.end(); ++I)
   {
      Result(j++, *I) = 1.0;
   }
   return Result;
}

SimpleOperator
construct_transfer_matrix(HermitianProxy<GenericMPO> const& A, GenericMPO const& B)
{
   CHECK_EQUAL(A.base().size(), B.size());

   SimpleOperator Result = local_inner_tensor_prod(herm(A.base()[0]), B[0]);
   for (unsigned i = 1; i < B.size(); ++i)
   {
      Result = Result *  local_inner_tensor_prod(herm(A.base()[i]), B[i]);
   }
   return Result;
}

// classification

OperatorClassification::OperatorClassification()
   : Factor_(0.0), Product_(false), Unitary_(false),
     Identity_(false), PropUnitary_(false), PropIdentity_(false), Null_(false)
{
}

bool OperatorClassification::is_null() const
{
   return Null_;
}

bool OperatorClassification::is_product() const
{
   return Product_;
}

bool OperatorClassification::is_unitary() const
{
   return Unitary_;
}

bool OperatorClassification::is_prop_unitary() const
{
   return PropUnitary_;
}

bool OperatorClassification::is_prop_identity() const
{
   return PropIdentity_;
}

bool OperatorClassification::is_identity() const
{
   return Identity_;
}

bool OperatorClassification::is_unclassified() const
{
   return !Product_ && !Null_;
}

std::complex<double> OperatorClassification::factor() const
{
   return Factor_;
}

std::ostream& operator<<(std::ostream& out, OperatorClassification const& Class)
{
   out << "null: " << Class.is_null() << '\n';
   out << "product: " << Class.is_product() << '\n';
   out << "unitary: " << Class.is_unitary() << '\n';
   out << "prop_unitary: " << Class.is_prop_unitary() << '\n';
   out << "prop_identity: " << Class.is_prop_identity() << '\n';
   out << "complex_identity: " << Class.is_complex_identity() << '\n';
   out << "identity: " << Class.is_identity() << '\n';
   out << "factor: " << Class.factor() << '\n';
   return out;
}

OperatorClassification classify(OperatorComponent c, double UnityEpsilon)
{
   OperatorClassification Result;

   // Early return if the component is null
   if (c.is_null())
   {
      Result.Null_ = true;
      return Result;
   }

   bool IsPropIdentity = true;  // true if the operator is proportional to identity
   bool IsPropUnitary = true;   // true if the operator is proportional to a unitary operator
   std::complex<double> Factor  = 1.0;

   // firstly, check to see if it is 1x1
   if (c.Basis1().size() != 1 || c.Basis2().size() != 1)
      return Result;  // default constructed return is unclassified

   if (IsPropUnitary)
   {
      SimpleRedOperator X = c(0,0);

      if (IsPropIdentity)
      {
         if (X.Basis1() != X.Basis2() || !is_pure_scalar(X))
            IsPropIdentity = false;
         else
         {
            std::complex<double> x = PropIdent(X.scalar(), UnityEpsilon);
            if (x == 0.0)
               IsPropIdentity = false;
            else
               Factor *= x;
         }
      }

      if (!IsPropIdentity)
      {
         // is it unitary?
         std::complex<double> x = PropIdent(scalar_prod(X, herm(X)), UnityEpsilon);
         std::complex<double> y = PropIdent(scalar_prod(herm(X), X), UnityEpsilon);

         if (x == 0.0 || y == 0.0)
         {
            IsPropUnitary = false;
         }
         else
         {
            Factor *= std::sqrt(x);
         }
      }
   }

   Result.Product_ = true;
   if (IsPropUnitary)
   {
      Result.PropUnitary_ = true;

      if (IsPropIdentity)
      {
         Result.PropIdentity_ = true;
         Result.Identity_ = LinearAlgebra::norm_frob_sq(Factor - std::complex<double>(1.0, 0))
            < UnityEpsilon*UnityEpsilon;

         // if we claim to be an identity operator, we might as well make it exact
         if (Result.Identity_)
            Factor = 1.0;

         Result.Unitary_ = LinearAlgebra::norm_frob_sq(norm_frob(Factor) - 1.0) < UnityEpsilon*UnityEpsilon;
         Result.Factor_ = Factor;
      }
      else
      {
         Result.Unitary_ = LinearAlgebra::norm_frob_sq(norm_frob(Factor) - 1.0) < UnityEpsilon*UnityEpsilon;
         Result.Factor_ = Factor;
      }
   }
   return Result;
}

OperatorClassification classify(GenericMPO const& Op, double UnityEpsilon)
{
   OperatorClassification Result;

   // Early return if the operator is null
   if (Op.is_null())
   {
      Result.Null_ = true;
      return Result;
   }

   bool IsPropIdentity = true;  // true if the operator is proportional to identity
   bool IsPropUnitary = true;   // true if the operator is proportional to a unitary operator
   std::complex<double> Factor  = 1.0;

   for (unsigned i = 0; i < Op.size(); ++i)
   {
      // firstly, check to see if it is 1x1
      if (Op[i].Basis1().size() != 1 || Op[i].Basis2().size() != 1)
         return Result;  // default constructed return is unclassified

      if (IsPropUnitary)
      {
         SimpleRedOperator X = Op[i](0,0);

         if (IsPropIdentity)
         {
            if (X.Basis1() != X.Basis2() || !is_pure_scalar(X))
               IsPropIdentity = false;
            else
            {
               std::complex<double> x = PropIdent(X.scalar(), UnityEpsilon);
               if (x == 0.0)
                  IsPropIdentity = false;
               else
                  Factor *= x;
            }
         }

         if (!IsPropIdentity)
         {
            // is it unitary?
            std::complex<double> x = PropIdent(scalar_prod(X, herm(X)), UnityEpsilon);
            std::complex<double> y = PropIdent(scalar_prod(herm(X), X), UnityEpsilon);

            if (x == 0.0 || y == 0.0)
            {
               IsPropUnitary = false;
            }
            else
            {
               Factor *= std::sqrt(x);
            }
         }

      }
   }

   Result.Product_ = true;
   if (IsPropUnitary)
   {
      Result.PropUnitary_ = true;

      if (IsPropIdentity)
      {
         Result.PropIdentity_ = true;
         Result.Identity_ = LinearAlgebra::norm_frob_sq(Factor - std::complex<double>(1.0, 0))
            < UnityEpsilon*UnityEpsilon;

         // if we claim to be an identity operator, we might as well make it exact
         if (Result.Identity_)
            Factor = 1.0;

         Result.Unitary_ = LinearAlgebra::norm_frob_sq(norm_frob(Factor) - 1.0) < UnityEpsilon*UnityEpsilon;
         Result.Factor_ = Factor;
      }
      else
      {
         Result.Unitary_ = LinearAlgebra::norm_frob_sq(norm_frob(Factor) - 1.0) < UnityEpsilon*UnityEpsilon;
         Result.Factor_ = Factor;
      }
   }

   return Result;
}

std::vector<BasisList>
ExtractLocalBasis1(GenericMPO const& Op)
{
   std::vector<BasisList> Result(Op.size());
   for (unsigned i = 0; i < Op.size(); ++i)
   {
      Result[i] = Op[i].LocalBasis1();
   }
   return Result;
}

std::vector<BasisList>
ExtractLocalBasis2(GenericMPO const& Op)
{
   std::vector<BasisList> Result(Op.size());
   for (unsigned i = 0; i < Op.size(); ++i)
   {
      Result[i] = Op[i].LocalBasis2();
   }
   return Result;
}

void optimize(GenericMPO& Op)
{
   if (Op.size() < 2)
      return;

   bool Reduced = true; // flag to indicate that we reduced a dimension
   // loop until we do a complete sweep with no reduction in dimensions
   while (Reduced)
   {
      Reduced = false;

      // Working left to right, optimize the Basis2
      SimpleOperator T = TruncateBasis2(Op.front());
      if (T.size1() != T.size2())
         Reduced = true;
      for (unsigned i = 1; i < Op.size()-1; ++i)
      {
         Op[i] = T * Op[i];
         T = TruncateBasis2(Op[i]);
         if (T.size1() != T.size2())
            Reduced = true;
      }
      Op.back() = T * Op.back();

      // Working right to left, optimize Basis1
      T = TruncateBasis1(Op.back());
      if (T.size1() != T.size2())
         Reduced = true;
      for (int i = Op.size()-2; i >= 1; --i)
      {
         Op[i] = Op[i] * T;
         T = TruncateBasis1(Op[i]);
         if (T.size1() != T.size2())
            Reduced = true;
      }
      Op.front() = Op.front() * T;
   }
}

BasisList ReverseBasis(BasisList const& B)
{
   std::vector<QuantumNumbers::QuantumNumber> BasisRev(B.begin(), B.end());
   std::reverse(BasisRev.begin(), BasisRev.end());
   std::transform(BasisRev.begin(), BasisRev.end(), BasisRev.begin(), [](auto x) { return adjoint(x); });
   return BasisList(B.GetSymmetryList(), BasisRev.begin(), BasisRev.end());
}

GenericMPO gauge_flip(GenericMPO const& Op)
{
   GenericMPO Result(Op.size());
   BasisList B1 = ReverseBasis(Op.Basis1());
   int i = 0;
   for (auto const& x : Op)
   {
      BasisList B2 = ReverseBasis(x.Basis2());
      OperatorComponent RevOp(x.LocalBasis1(), x.LocalBasis2(), B2, B1);
      for (OperatorComponent::const_iterator I = iterate(x); I; ++I)
      {
         for (OperatorComponent::const_inner_iterator J = iterate(I); J; ++J)
         {
            RevOp(B2.size()-J.index2()-1, B1.size()-J.index1()-1) = *J;
         }
      }
      Result[i] = RevOp;
      B1 = B2;
      ++i;
   }
   Result.debug_check_structure();
   return Result;
}
