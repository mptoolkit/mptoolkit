// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mpo/generic_mpo.h
//
// Copyright (C) 2012-2022 Ian McCulloch <ian@qusim.net>
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

// GenericMPO represents an MPO that is in no particular form
// It contains a LatticeCommute, so that we can determine the
// Jordan-Wigner string that should be used when extending the
// operator onto a lattice.

#if !defined(MPTOOLKIT_MPO_GENERIC_MPO_H)
#define MPTOOLKIT_MPO_GENERIC_MPO_H

#include "operator_component.h"
#include <vector>

class GenericMPO
{
   private:
      typedef std::vector<OperatorComponent> DataType;

   public:
      typedef OperatorComponent value_type;
      typedef DataType::iterator iterator;
      typedef DataType::const_iterator const_iterator;
      typedef value_type::basis1_type basis1_type;
      typedef value_type::basis2_type basis2_type;

      GenericMPO() {}

      explicit GenericMPO(int Size) : Data_(Size) {}

      explicit GenericMPO(OperatorComponent const& x) : Data_(1, x) {}

      // Size repeated copies of x
      GenericMPO(int Size, OperatorComponent const& x) : Data_(Size, x) {}

      // from an iterator
      template <typename InIter>
      GenericMPO(InIter first, InIter last) : Data_(first, last) {}

      bool empty() const { return Data_.empty(); }
      std::size_t size() const { return Data_.size(); }

      bool is_null() const;

      iterator begin() { return Data_.begin(); }
      iterator end() { return Data_.end(); }

      const_iterator begin() const { return Data_.begin(); }
      const_iterator end() const { return Data_.end(); }

      value_type& operator[](int x) { return Data_[x]; }
      value_type const& operator[](int x) const { return Data_[x]; }

      value_type const& front() const { return Data_.front(); }
      value_type& front() { return Data_.front(); }

      value_type const& back() const { return Data_.back(); }
      value_type& back() { return Data_.back(); }

      basis1_type Basis1() const { DEBUG_CHECK(!this->empty()); return Data_.front().Basis1(); }
      basis2_type Basis2() const { DEBUG_CHECK(!this->empty()); return Data_.back().Basis2(); }

      QuantumNumbers::SymmetryList GetSymmetryList() const { return Data_[0].GetSymmetryList(); }

      // Return the local basis at the n'th site
      BasisList const& LocalBasis1(int n) const
      { return Data_[n].LocalBasis1(); }
      BasisList const& LocalBasis2(int n) const
      { return Data_[n].LocalBasis2(); }

      // returns the list of local hilbert spaces for this operator
      std::vector<BasisList> LocalBasis1List() const;
      std::vector<BasisList> LocalBasis2List() const;

      std::vector<OperatorComponent> const& data() const { return Data_; }

      void check_structure() const;

      void debug_check_structure() const;

   private:
      DataType Data_;

   friend PStream::opstream& operator<<(PStream::opstream& out, GenericMPO const& op);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, GenericMPO& op);
};

std::ostream&
operator<<(std::ostream& out, GenericMPO const& op);

PStream::opstream&
operator<<(PStream::opstream& out, GenericMPO const& op);

PStream::ipstream&
operator>>(PStream::ipstream& in, GenericMPO& op);

// remove unused matrix elements
void cull_unused_elements(GenericMPO& Op);

// As an alternative to cull_unused_elements(), we can use a mask vector which indicates which components are unused.
// This is a 2D array, indexed by the bond number of the MPO, and then the index into the auxiliary basis.
// initialize_mask sets
typedef std::vector<std::vector<int> > MPOMaskType;
void initialize_mask(GenericMPO const& Op, MPOMaskType& Mask);

// Updates a mask to exclude matrix elements that won't be used because they would be zero or don't contribute
// to the result.
void mask_unused_elements(GenericMPO const& Op, std::vector<std::vector<int> >& Mask);

// Construct a mask that singles out a particular column of the MPO
MPOMaskType mask_column(GenericMPO const& Op, int Col);

// Construct a mask that singles out a particular row of the MPO
MPOMaskType mask_row(GenericMPO const& Op, int Row);

// Does a N-1 coarse graining of an operator.  The length must be a multiple of N
GenericMPO coarse_grain(GenericMPO const& Op, int N);

// Coarse-grains a section of an MPO into a single site.
GenericMPO coarse_grain_range(GenericMPO const& Op, int beg, int end);

// constructs the transfer operator as
// prod_i local_inner_tensor_prod(herm(A.base()[i]), B[i])
SimpleOperator
construct_transfer_matrix(HermitianProxy<GenericMPO> const& A, GenericMPO const& B);

struct OperatorClassification
{
   // indicates that the operator is zero
   bool is_null() const;

   // indicates that the operator a product state, ie a product of 1x1 MPO's
   bool is_product() const;

   // indicates that the operator is a unitary product state, ie a string operator
   bool is_unitary() const;

   // indicates that the operator is proportional to a unitary product state,
   // up to some factor
   bool is_prop_unitary() const;

   // indicates that the operator is proportional to the identity operator
   bool is_prop_identity() const;

   // returns true if the operator is the identity multiplied by a complex phase factor of magnitude 1
   bool is_complex_identity() const
   {
      return this->is_prop_identity() && norm_frob(norm_frob(this->factor())-1.0) < 1E-12;
   }

   // indicates that the operator is equal to the identity
   bool is_identity() const;

   // returns true only if the operator fits into no other classification
   bool is_unclassified() const;

   // for operators that are proportional to the identity, returns the factor
   std::complex<double> factor() const;

   // private use only
   std::complex<double> Factor_;
   bool Product_;
   bool Unitary_;
   bool Identity_;
   bool PropUnitary_;
   bool PropIdentity_;
   bool Null_;

   OperatorClassification();
};

std::ostream& operator<<(std::ostream& out, OperatorClassification const& Class);

OperatorClassification classify(GenericMPO const& Op, double UnityEpsilon);

inline
OperatorClassification classify(GenericMPO const& Op)
{
   return classify(Op, DefaultClassifyUnityEpsilon);
}

OperatorClassification classify(OperatorComponent c, double UnityEpsilon);

inline
OperatorClassification classify(OperatorComponent c)
{
   return classify(c, DefaultClassifyUnityEpsilon);
}

// plus various functions for acting on states etc

// TODO: find a better location for this function
// Construct an operator that projects onto a given subset of a basis.
SimpleOperator make_projector_onto(BasisList const& Basis, std::set<int> const& Onto);

// extract the local basis at each site of the MPO
std::vector<BasisList>
ExtractLocalBasis1(GenericMPO const& Op);

std::vector<BasisList>
ExtractLocalBasis2(GenericMPO const& Op);

// Reverse the direction of the  rows and coumns of an MPO.  SO
// (A B)
// (0 C)
// turns into
// (C B)
// (- A)
GenericMPO gauge_flip(GenericMPO const& Op);

namespace LinearAlgebra
{

template <>
struct interface<GenericMPO>
{
   typedef void type;
};

template <>
struct Herm<GenericMPO>
{
   typedef GenericMPO const& argument_type;
   typedef HermitianProxy<GenericMPO> result_type;

   result_type operator()(argument_type x) const
   { return result_type(x); }
};

template <>
struct Conj<GenericMPO>
{
   typedef GenericMPO const& argument_type;
   typedef GenericMPO result_type;

   result_type operator()(argument_type x) const
   {
      GenericMPO Result(x);
      for (GenericMPO::iterator I = Result.begin(); I != Result.end(); ++I)
      {
         *I = conj(*I);
      }
      return Result;
   }
};

template <>
struct Adjoint<GenericMPO>
{
   typedef GenericMPO const& argument_type;
   typedef GenericMPO result_type;

   result_type operator()(argument_type x) const
   {
      GenericMPO Result(x);
      for (GenericMPO::iterator I = Result.begin(); I != Result.end(); ++I)
      {
         *I = adjoint(*I);
      }
      return Result;
   }
};

template <>
struct InvAdjoint<GenericMPO>
{
   typedef GenericMPO const& argument_type;
   typedef GenericMPO result_type;

   result_type operator()(argument_type x) const
   {
      GenericMPO Result(x);
      for (GenericMPO::iterator I = Result.begin(); I != Result.end(); ++I)
      {
         *I = inv_adjoint(*I);
      }
      return Result;
   }
};

} // namespace LinearAlgebra

#include "generic_mpo.cc"

#endif
