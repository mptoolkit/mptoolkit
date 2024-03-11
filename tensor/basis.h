// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// tensor/basis.h
//
// Copyright (C) 2004-2024 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2023 Jesse Osborne <j.osborne@uqconnect.edu.au>
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
/* -*- C++ -*- $Id$

  Extends BasisList to keep track of the dimension of each subspace, hence
  the name 'VectorBasis', as each basis label represents a vector of states.

*/

#if !defined(MPTOOLKIT_TENSOR_BASIS_H)
#define MPTOOLKIT_TENSOR_BASIS_H

#include "quantumnumbers/symmetrylist.h"
#include "quantumnumbers/quantumnumber.h"
#include "linearalgebra/vector.h"
#include "linearalgebra/pstreamio.h"
#include "linearalgebra/sparsematrix.h"
#include "linearalgebra/diagonalmatrix.h"

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
      typedef QuantumNumber value_type;

      BasisList() {}

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

      typedef QuantumNumbers::QuantumNumberList::const_iterator const_iterator;

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

   friend BasisList delta_shift(BasisList const& b, QuantumNumber const& q);

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

std::ostream&
operator<<(std::ostream& out, std::vector<BasisList> const& u);

// Returns true if two vectors of BasisList are compatible.
// A compatible site list is one where we firstly make each vector
// the same size by repeating the unit cell as many times as necessary,
// and the lists are compatible if the individual BasisList's match exactly
// at each site.
bool
is_compatible(std::vector<BasisList> const& a, std::vector<BasisList> const& b);

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

// Get a map of the dimension of the basis per quantum number sector
inline
std::map<QuantumNumbers::QuantumNumber, int>
DimensionPerSector(BasisList const& b)
{
   std::map<QuantumNumbers::QuantumNumber, int> Result;
   for (int i = 0; i < b.size(); ++i)
      ++Result[b[i]];
   return Result;
}

// construct a mapping from the components with a particular quantum number,
// onto successive integers
std::map<int, int>
LinearizeQuantumNumberSubspace(BasisList const& b, QuantumNumbers::QuantumNumber const& q);

// Apply a shift operator to the basis.  This will fail if any of the shifts
// are not possible (eg, if it shifts beyond the highest weight rep)
BasisList delta_shift(BasisList const& Orig, QuantumNumbers::Projection const& p);

BasisList delta_shift(BasisList const& b, QuantumNumber const& q);


//
// VectorBasis
//
// Extends the concept of BasisList to maintain a dimension for each subspace.
//

class VectorBasis
{
   public:
      typedef BasisList::value_type value_type;

      VectorBasis() {}

      explicit VectorBasis(SymmetryList const& sl);
      explicit VectorBasis(BasisList const& Basis);

      // Constructor from a BasisList and a container of dimensions
      template <typename FwdIter>
      VectorBasis(BasisList const& Basis, FwdIter first, FwdIter last);

      // Constructor from a container of std::pair<QuantumNumber, integer>
      template <typename FwdIter>
      VectorBasis(FwdIter first, FwdIter last);

      // Version of the above which sets the symmetry list: this allows the container to be empty.
      template <typename FwdIter>
      VectorBasis(QuantumNumbers::SymmetryList const& S, FwdIter first, FwdIter last);

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
      VectorBasis(BasisList const& b, LinearAlgebra::Vector<int> const& dim)
         : Basis_(b), Dimension_(dim) {}

      BasisList Basis_;
      LinearAlgebra::Vector<int> Dimension_;

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

inline
std::map<QuantumNumbers::QuantumNumber, int>
DimensionPerSector(VectorBasis const& b)
{
   std::map<QuantumNumbers::QuantumNumber, int> Result;
   for (int i = 0; i < b.size(); ++i)
      Result[b[i]] += b.dim(i);
   return Result;
}

// Given some opearator over basis (B1,B2), return the rank of each quantum number sector, as the
// map from each quantum number -> minimum of the dimension of that sector in B1 and B2.
// There is also a version that uses a TensorOperator, defined in tensor.h, that also looks at the
// structure of the tensor.
inline
std::map<QuantumNumbers::QuantumNumber, int>
RankPerSector(VectorBasis const& B1, VectorBasis const& B2)
{
   std::map<QuantumNumbers::QuantumNumber, int> Basis1Size, Basis2Size;
   for (int i = 0; i < B1.size(); ++i)
      Basis1Size[B1[i]] += B1.dim(i);
   for (int j = 0; j < B2.size(); ++j)
      Basis2Size[B2[j]] += B2.dim(j);
   std::map<QuantumNumbers::QuantumNumber, int> Result;
   for (auto const& q : Basis1Size)
   {
      if (q.second > 0 && Basis2Size[q.first] > 0)
         Result[q.first] = std::min(q.second, Basis2Size[q.first]);
   }
   return Result;
}

// Returns the dimension of the quantum number sectors in B1, but with sectors that don't exist in B2
// removed.
inline
std::map<QuantumNumbers::QuantumNumber, int>
CullMissingSectors(VectorBasis const& B1, VectorBasis const& B2)
{
   std::set<QuantumNumbers::QuantumNumber> QuantumNumbersInB2 = QuantumNumbersInBasis(B2);
   std::map<QuantumNumbers::QuantumNumber, int> Result;
   for (int i = 0; i < B1.size(); ++i)
   {
      if (QuantumNumbersInB2.count(B1[i]) > 0)
         Result[B1[i]] += B1.dim(i);
   }
   return Result;
}

// Given two vector bases B1 and B2, construct a third basis that for each quantum number q that is in
// B1 and B2, has dimension equal to min(d1,d2), where d1,d2 are the dimensions of the q sector in B1 and B2
// respectively.  The second and third components of the tuple are arrays of size B1.size() and B2.size() respectively,
// which contain the index into the new basis of the corresponding sector of B1 and B2, or -1 if the new basis
// doesn't contain this sector.
// B1 abd B2 must be regular.
std::tuple<VectorBasis, std::vector<int>, std::vector<int>>
MakeMinimalBasis(VectorBasis const& B1, VectorBasis const& B2);

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

template <typename T, typename B1, typename B2>
struct MakeZeroImpl {};

// NOTE: T is not deduced here
template <typename T, typename B1, typename B2>
typename MakeZeroImpl<T, B1, B2>::result_type
make_zero(B1 const& b1, B2 const& b2, int i, int j)
{
   return MakeZeroImpl<T, B1, B2>(b1, b2, i, j);
}

template <typename T>
struct MakeZeroImpl<T, VectorBasis, VectorBasis>
{
   typedef T result_type;
   T operator()(VectorBasis const& b1, VectorBasis const& b2, int i, int j) const
   {
      return LinearAlgebra::SparseMatrix<double>(b1.dim(i), b2.dim(j));
   }
};

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
      //      return LinearAlgebra::identity_matrix<double>(b.dim(i));
      return LinearAlgebra::DiagonalMatrix<double>(b.dim(i), b.dim(i), 1.0);
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
