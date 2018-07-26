// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// tensor/basis.h
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

#if !defined(MPTOOLKIT_TENSOR_BASIS_H)
#define MPTOOLKIT_TENSOR_BASIS_H

#include "quantumnumbers/symmetrylist.h"
#include "quantumnumbers/quantumnumber.h"

namespace Tensor
{

using QuantumNumbers::SymmetryList;
using QuantumNumbers::QuantumNumber;
using QuantumNumbers::QuantumNumberList;
using QuantumNumbers::Projection;

//
// BasisList
//
// This holds the quantum numbers of each index of a tensor operator.
//

class BasisList
{
   public:
      using value_type     = QuantumNumber;
      using const_iterator = QuantumNumbers::QuantumNumberList::const_iterator;

      BasisList() noexcept {}

      ~BasisList() noexcept = default;

      explicit BasisList(QuantumNumbers::SymmetryList const& S) : S_(S) {}

      // construction of a singleton BasisList from a quantum number
      explicit BasisList(QuantumNumbers::QuantumNumber const& q) : S_(q.GetSymmetryList()), Q_(1,q) {}

      // Construction from a list of quantum numbers.  The list MUST be non-empty
      // so we can set the SymmetryList
      template <typename FwdIter>
      BasisList(FwdIter first, FwdIter last);

      // Construction from a list of quantum numbers, and set the symmetry list.
      // This allows the list to be empty.
      template <typename FwdIter>
      BasisList(QuantumNumbers::SymmetryList const& S, FwdIter first, FwdIter last);

      BasisList(BasisList const& Other) = default;
      BasisList(BasisList&& Other) noexcept = default;

      BasisList& operator=(BasisList const& Other) = default;
      BasisList& operator=(BasisList&& Other) noexcept = default;

      const_iterator begin() const { return Q_.begin(); }
      const_iterator end() const { return Q_.end(); }

      value_type const& operator[](int x) const { return Q_[x]; }

      // returns the first quantum number in the basis
      value_type const& front() const { return Q_.front(); }

      // returns the last quantum number in the basis
      value_type const& back() const { return Q_.back(); }

      std::size_t size() const { return Q_.size(); }

      int dim(int s) const { DEBUG_CHECK(s >= 0 && s < int(this->size())); return 1; }

      bool is_null() const { return S_.is_null(); }

      // returns true if this basis is empty
      bool is_empty() const { return this->size() == 0; }

      void push_back(QuantumNumber const& q)
         { DEBUG_PRECONDITION_EQUAL(q.GetSymmetryList(), S_); Q_.push_back(q); }

      // returns true if this BasisList contains only one element, which transforms as the scalar quantum number
      bool is_identity() const;

      // returns true if this BasisList is 'regular'; that is, each quantum number occurs at most once in the basis
      bool is_regular() const;

      int total_degree() const;

      // finds the first occurance of quantum number q in the basis, or -1 if no such q
      int find_first(QuantumNumber const& q) const;

      // Finds the next occurance of quantum number q in the basis beyond n, or -1 if no such q
      int find_next(QuantumNumber const& q, int n) const;

      SymmetryList const& GetSymmetryList() const { return S_; }

      void CoerceSymmetryList(SymmetryList const& sl);

      void delta_shift(QuantumNumber const& q);

   private:
      QuantumNumbers::SymmetryList S_;
      QuantumNumbers::QuantumNumberList Q_;

   friend bool operator==(BasisList const& b1, BasisList const& b2)
      { return b1.Q_ == b2.Q_; }

   friend bool operator!=(BasisList const& b1, BasisList const& b2)
      { return b1.Q_ != b2.Q_; }

   friend PStream::opstream& operator<<(PStream::opstream& out, BasisList const& B);
      friend PStream::ipstream& operator>>(PStream::ipstream& in, BasisList& B);
   friend BasisList CoerceSymmetryList(BasisList const& b, SymmetryList const& sl)
   __attribute__((warn_unused_result));

   friend void CoerceSymmetryListInPlace(BasisList& b, SymmetryList const& sl);
};

typedef BasisList SimpleBasis;  // for backwards compatibility

// Helper to construct a basis with a single state, the identity quantum number
BasisList
make_vacuum_basis(SymmetryList const& S);

// Helper to construct a basis with a single state of the specified quantum number
BasisList
make_single_basis(QuantumNumbers::QuantumNumber const& q);


std::ostream& operator<<(std::ostream& out, BasisList const& b);

BasisList adjoint(BasisList const& b);

std::string show_projections(BasisList const& B);


inline
void
BasisList::CoerceSymmetryList(SymmetryList const& sl)
{
   S_ = sl;
   Q_.CoerceSymmetryList(sl);
}

inline
BasisList
CoerceSymmetryList(BasisList const& b, SymmetryList const& sl)
{
   BasisList Result;
   Result.S_ = sl;
   Result.Q_ = CoerceSymmetryList(b.Q_, sl);
   return Result;
}

inline
void
CoerceSymmetryListInPlace(BasisList& b, SymmetryList const& sl)
{
   b.CoerceSymmetryList(sl);
}

//BasisList RenameSymmetry(BasisList const& BL, SymmetryList const& NewSL);

inline
std::set<QuantumNumbers::QuantumNumber>
QuantumNumbersInBasis(BasisList const& b)
{
   return std::set<QuantumNumbers::QuantumNumber>(b.begin(), b.end());
}

// construct a mapping from the components with a particular quantum number,
// onto successive integers
std::map<int, int>
LinearizeQuantumNumberSubspace(BasisList const& b, QuantumNumbers::QuantumNumber const& q);

// Apply a shift operator to the basis.  This will fail if any of the shifts
// are not possible (eg, if it shifts beyond the highest weight rep)
BasisList delta_shift(BasisList const& Orig, QuantumNumbers::Projection const& p);

//
// VectorBasis
//
// Extends the concept of BasisList to maintain a dimension for each subspace.
//

class VectorBasis
{
   public:
      typedef BasisList::value_type value_type;

      VectorBasis() noexcept {}

      static_assert(std::is_nothrow_move_constructible<BasisList>::value, "");
      static_assert(std::is_nothrow_move_constructible<std::vector<int> >::value, "");


      explicit VectorBasis(SymmetryList const& sl);
      explicit VectorBasis(BasisList Basis);

      // Constructor from a BasisList and a container of dimensions
      template <typename FwdIter>
      VectorBasis(BasisList const& Basis, FwdIter first, FwdIter last);

      // Constructor from a container of std::pair<QuantumNumber, integer>
      template <typename FwdIter>
      VectorBasis(FwdIter first, FwdIter last);

      VectorBasis(VectorBasis const& Other) = default;
      VectorBasis(VectorBasis&& Other) noexcept = default;
      VectorBasis& operator=(VectorBasis const& Other) = default;
      VectorBasis& operator=(VectorBasis&& Other) noexcept = default;

      void push_back(QuantumNumber const& q, int Dimension);

      void push_back(QuantumNumber const& q)
      { this->push_back(q, 0); }

      // returns the number of subspaces in the basis
      std::size_t size() const { return Basis_.size(); }

      bool is_null() const { return Basis_.is_null(); }

      // returns true if this basis is empty
      bool is_empty() const { return Basis_.size() == 0; }

      // returns true if the basis is a 'vacuum basis', that is,
      // a 1-dimensional basis in the scalar sector
      bool is_vacuum() const { return Basis_.size() == 1 && is_scalar(Basis_[0]) &&
	    Dimension_[0] == 1; }

      value_type const& operator[](int s) const { return Basis_[s]; }

      void set_dim(int s, int d) { Dimension_[s] = d; }

      int dim(int s) const { return Dimension_[s]; }

      // returns the total number of states in the basis.
      int total_dimension() const;

      // returns the degree of the group representation, summed over every basis state.
      int total_degree() const;

      SymmetryList const& GetSymmetryList() const { return Basis_.GetSymmetryList(); }

      // finds the first occurance of quantum number q in the basis, or -1 if no such q
      int find_first(QuantumNumber const& q) const
      { return Basis_.find_first(q); }

      // Finds the next occurance of quantum number q in the basis beyond n, or -1 if no such q
      int find_next(QuantumNumber const& q, int n) const
      { return Basis_.find_next(q,n); }

      BasisList const& Basis() const { return Basis_; }
      BasisList& Basis() { return Basis_; }

      void CoerceSymmetryList(SymmetryList const& SL);

      // implicit conversion to BasisList **OUCH** why is this here?
      operator BasisList const&() const { return Basis_; }

      void delta_shift(QuantumNumber const& q);

 private:
      VectorBasis(BasisList b, std::vector<int> dim)
         : Basis_(std::move(b)), Dimension_(std::move(dim)) {}

      BasisList Basis_;
      std::vector<int> Dimension_;

   friend PStream::opstream& operator<<(PStream::opstream& out, VectorBasis const& Bi);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, VectorBasis& Bi);

   friend bool operator==(VectorBasis const& x, VectorBasis const& y)
      { return x.Basis_ == y.Basis_ && x.Dimension_ == y.Dimension_; }
   friend bool operator!=(VectorBasis const& x, VectorBasis const& y)
      { return !operator==(x,y); }

   friend VectorBasis CoerceSymmetryList(VectorBasis const& b, SymmetryList const& sl)
   __attribute__((warn_unused_result));


   friend void CoerceSymmetryListInPlace(VectorBasis& b, SymmetryList const& sl);

   friend VectorBasis delta_shift(VectorBasis const& Orig,
                                 QuantumNumbers::Projection const& p);

   friend VectorBasis delta_shift(VectorBasis const& Orig,
                                 QuantumNumbers::QuantumNumber const& q);

   friend VectorBasis RenameSymmetry(VectorBasis const& BL, SymmetryList const& NewSL);

   friend VectorBasis adjoint(VectorBasis const& b);
};

inline
void
VectorBasis::CoerceSymmetryList(SymmetryList const& sl)
{
   Basis_.CoerceSymmetryList(sl);
}

inline
void
CoerceSymmetryListInPlace(VectorBasis& b, SymmetryList const& sl)
{
   b.CoerceSymmetryList(sl);
}

std::ostream& operator<<(std::ostream& out, VectorBasis const& B);

VectorBasis adjoint(VectorBasis const& b);

std::string show_projections(VectorBasis const& B);

inline
std::set<QuantumNumbers::QuantumNumber>
QuantumNumbersInBasis(VectorBasis const& b)
{
   return std::set<QuantumNumbers::QuantumNumber>(b.Basis().begin(), b.Basis().end());
}

VectorBasis RenameSymmetry(VectorBasis const& BL, SymmetryList const& NewSL);

// Apply a shift operator to the basis.  This will fail if any of the shifts
// are not possible (eg, if it shifts beyond the highest weight rep)
VectorBasis delta_shift(VectorBasis const& Orig, QuantumNumbers::Projection const& p);

#if 0
// Apply a shift operation, where q is a degree 1 rep
inline
VectorBasis
delta_shift(VectorBasis const& Orig, QuantumNumbers::QuantumNumber const& q)
{
   QuantumNumbers::ProjectionList PL = enumerate_projections(q);
   DEBUG_PRECONDITION_EQUAL(PL.size(), 1);
   return delta_shift(Orig, PL[0]);
}

inline
VectorBasis
delta_shift(VectorBasis const& Orig, QuantumNumbers::QuantumNumber const& q)
{
   return delta_shift(Orig, q);
}
#endif

//
// make_zero
//

//
// make_identity
//

template <typename T, typename B>
struct MakeIdentityImpl {};

// NOTE: T is not deduced here
template <typename T, typename B>
typename MakeIdentityImpl<T, B>::result_type
make_identity(B const& b, std::size_t i)
{
   return MakeIdentityImpl<T, B>()(b, i);
}

template <typename T>
struct MakeIdentityImpl<T, VectorBasis>
{
   typedef T result_type;
   T operator()(VectorBasis const& b, int i) const
   {
      return T::make_identity(b.dim(i));
   }
};

template <typename T>
struct MakeIdentityImpl<T, BasisList>
{
   typedef T result_type;
   T operator()(BasisList const& b, int i) const
   {
      return T(1.0);
   }
};

} // namespace Tensor

#include "basis.cc"

#endif
