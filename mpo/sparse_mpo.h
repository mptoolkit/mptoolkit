// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mpo/sparse_mpo.h
//
// Copyright (C) 2013-2020 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

// SparseMPO represents a generic MPO which is stored
// in sparse form, where only elements that are not diagonal/identity
// are actually stored.

#if !defined(MPTOOLKIT_MPO_COMPRESSED_MPO_H)
#define MPTOOLKIT_MPO_COMPRESSED_MPO_H

#include "generic_mpo.h"
#include <map>

class SparseMPO
{
   public:
      // The storage is a map of site indices into the elements of a GenericMPO
      using MapType = std::map<int, int>;

      using value_type = OperatorComponent;

      class const_iterator;

      using base_iterator       = DataType::iterator;
      using const_base_iterator = DataType::const_iterator;

      SparseMPO() {}

      explicit SparseMPO(SiteBasisList const& SList);

      bool empty() const { return SiteList_.empty(); }

      std::size_t size() const { return SiteList_.size(); }

      bool is_null() const { return Data_.empty(); }

      const_iterator begin() const;
      const_iterator end() const;

      OperatorComponent operator[](int x) const;

      iterator base_begin() { return Data_.begin(); }
      iterator base_end() { return Data_.end(); }

      // returns the non-sparse part of the MPO
      GenericMPO const& values() const { return Data_; }
      GenericMPO& values() { return Data_; }

      BasisList Basis1() const;
      BasisList Basis2() const;

      SymmetryList GetSymmetryList() const { return SiteList_.GetSymmetryList(); }

      // Return the local basis at the n'th site
      BasisList const& LocalBasis1(int n) const;

      BasisList const& LocalBasis2(int n) const;

      // returns the list of local hilbert spaces for this operator
      SiteBasisList LocalBasis1List() const;
      SiteBasisList LocalBasis2List() const;

      std::vector<OperatorComponent> const& data() const { return Data_; }

      // iteration over the base GenericMPO
      const_iterator base_begin() const { return Data_.begin(); }
      const_iterator base_end() const { return Data_.end(); }

      // iteration over the sparse mapping
      MapType::const_iterator map_begin() const { return Map_.begin(); }
      MapType::const_iterator map_end() const { return Map_.end(); }

      void check_structure() const;

      void debug_check_structure() const;

   private:
      SiteBasisList SiteList_;
      MapType Map_;
      DataType Data_;

   friend PStream::opstream& operator<<(PStream::opstream& out, SparseMPO const& op);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, SparseMPO& op);
};

class SparseMPO::const_iterator
{
   public:
      const_iterator();

      using value_type = OperatorComponent;

      void operator++();

      void operator--();

      value_type const& operator*();

      value_type const* operator->();

   protected:
      bool is_at_end() const { return Loc == Container->size(); }
      bool is_at_begin() const { return Loc == 0; }

   private:
      const_iterator(SparseMPO const* Container_, const_base_iterator I_, MapType::const_iterator M_, SiteBasisList::const_iterator S_, int Loc_);

      SparseMPO const* Container;
      GenericMPO::const_iterator I;
      MapType::const_iterator M;
      SiteBasisList::const_iterator S;
      OperatorComponent C;
      int Loc;
}

SparseMPO::const_iterator::const_iterator() : Container(nullptr), Loc(0)
{
}

SparseMPO::const_iterator::const_iterator(SparseMPO const* Container_, GenericMPO::const_iterator I_, MapType::const_iterator M_, SiteBasisList::const_iterator S_, int Loc_)
   : Container(Container_), I(I_), M(M_), S(S_), Loc(Loc_)
{
}

void SparseMPO::const_iterator::operator++()
{
   C = OperatorComponent();
   ++Loc;
   ++S;
   if (M != Container->map_end() && Loc > M->first)
   {
      ++M;
      ++I;
   }
}

void SparseMPO::const_iterator::operator--()
{
   C = OperatorComponent();
   --Loc;
   --S;
   if (M != Container->map_end() && Loc < M->first && M != Container->map_begin())
   {
      --I;
      --M;
   }
}

OperatorComponent const&
SparseMPO::const_iterator::operator*()
{
   // if we're sitting on the site of a value, then return the value
   if (M != Container->map_end() && Loc == M->first)
      return *I;

   // if we've already initialized the operator component, then return it
   if (!C.is_null())
      return C;

   // otherwise we need to initialize C.  Is M pointing beyond the end of the last physical site?
   if (M == Container->map_end())
   {
      C = OperatorComponent::make_identity(*S, Container->Basis2());
      return C;
   }

   // at this point, M is a valid value which could be either to the left or right of Loc
   if (Loc < M->first)
   {
      C = OperatorComponent::make_identity(*S, I->Basis1())
   }
   else
   {
      DEBUG_CHECK(Loc > M->first);
      C = OperatorComponent::make_identity(*S, I->Basis2())
   }
   return C;
}

OperatorComponent const*
SparseMPO::const_iterator::operator->()
{
   return &this->operator*();
}

std::ostream&
operator<<(std::ostream& out, SparseMPO const& op);

PStream::opstream&
operator<<(PStream::opstream& out, SparseMPO const& op);

PStream::ipstream&
operator>>(PStream::ipstream& in, SparseMPO& op);

// remove unused matrix elements
void cull_unused_elements(SparseMPO& Op);

// Does a N-1 coarse graining of an operator.  The length must be a multiple of N
SparseMPO coarse_grain(SparseMPO const& Op, int N);

// Coarse-grains a section of an MPO into a single site.
SparseMPO coarse_grain_range(SparseMPO const& Op, int beg, int end);

// constructs the transfer operator as
// prod_i local_inner_tensor_prod(herm(A.base()[i]), B[i])
//SimpleOperator
//construct_transfer_matrix(HermitianProxy<SparseMPO> const& A, SparseMPO const& B);

// extract the local basis at each site of the MPO
std::vector<BasisList>
ExtractLocalBasis1(SparseMPO const& Op);

std::vector<BasisList>
ExtractLocalBasis2(SparseMPO const& Op);

namespace LinearAlgebra
{

template <>
struct interface<SparseMPO>
{
   typedef void type;
};

template <>
struct Herm<SparseMPO>
{
   typedef SparseMPO const& argument_type;
   typedef HermitianProxy<SparseMPO> result_type;

   result_type operator()(argument_type x) const
   { return result_type(x); }
};

template <>
struct Conj<SparseMPO>
{
   typedef SparseMPO const& argument_type;
   typedef SparseMPO result_type;

   result_type operator()(argument_type x) const
   {
      SparseMPO Result(x);
      for (SparseMPO::iterator I = Result.begin(); I != Result.end(); ++I)
      {
         *I = conj(*I);
      }
      return Result;
   }
};

template <>
struct Adjoint<SparseMPO>
{
   typedef SparseMPO const& argument_type;
   typedef SparseMPO result_type;

   result_type operator()(argument_type x) const
   {
      SparseMPO Result(x);
      for (SparseMPO::iterator I = Result.begin(); I != Result.end(); ++I)
      {
         *I = adjoint(*I);
      }
      return Result;
   }
};

template <>
struct InvAdjoint<SparseMPO>
{
   typedef SparseMPO const& argument_type;
   typedef SparseMPO result_type;

   result_type operator()(argument_type x) const
   {
      SparseMPO Result(x);
      for (SparseMPO::iterator I = Result.begin(); I != Result.end(); ++I)
      {
         *I = inv_adjoint(*I);
      }
      return Result;
   }
};

} // namespace LinearAlgebra

#endif
