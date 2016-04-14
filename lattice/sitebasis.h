/* -*- C++ -*- $Id$

  Version 2,
  Created 2000-06-16 Ian McCulloch

  Modified to use the sparse column format SparseMatrix rather than 
  the completely sparse std::map.  This is to partially merge the SiteOperator
  and BlockOperator concepts.  Ie, a SiteOperator is now essentially
  just the LinearSparsePart of a BlockOperator.
*/

#if !defined(SITEBASIS_H_F56734H8972NB78F3UY568)
#define SITEBASIS_H_F56734H8972NB78F3UY568

#include "tensor/tensor.h"
#include "tensor/tensorproduct.h"
#include "linearalgebra/matrix.h"
#include "common/numerics.h"
#include <map>
#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <boost/tuple/tuple.hpp>
#include <pheap/pvalueptr.h>

using QuantumNumbers::QuantumNumber;
using QuantumNumbers::Projection;
using QuantumNumbers::ProjectionList;
using QuantumNumbers::SymmetryList;
using LinearAlgebra::size_type;
using Tensor::BasisList;
using Tensor::ProductBasis;

//
// SiteBasis
//

class SiteBasis
{
   public:
      typedef std::pair<std::string, QuantumNumber> value_type;

      SiteBasis() : Label_(new LabelType()) {}

      // constructs an empty site basis
      SiteBasis(SymmetryList const& SList);

      SiteBasis(std::string const& SList);

      // Adds a quantum number to the basis.  Label is a unique string identifier,
      // and q is the quantum number of the state.  All states have a dimension 1,
      // thus the label uniquely identifies a single basis state.
      void push_back(std::string const& Label, QuantumNumber const& q);

      QuantumNumber const& qn(int i) const { return Basis_[i]; }

      value_type operator[](int i) const { return value_type((*Label_)[i], Basis_[i]); }

      size_type size() const { return Basis_.size(); }

      // returns the label corresponding to the given state.
      int Lookup(std::string const& Label) const;

      // returns the label corresponding to the state, or a negative number if no such state
      int LookupOrNeg(std::string const& Label) const;

      QuantumNumber const& qn(std::string const& Label) const
      { return this->qn(this->Lookup(Label)); }

      // Returns the string label of a given basis state
      std::string const& Label(int l) const { return (*Label_)[l]; }

      // implicit conversion to BasisList
      operator BasisList const&() const { return Basis_; }

      // explicit conversion, and backwards compatibility
      BasisList const& Basis() const { return Basis_; }

      SymmetryList const& GetSymmetryList() const { return Basis_.GetSymmetryList(); }

      bool operator==(SiteBasis const& S2) const
      { return Basis_ == S2.Basis_ && Label_ == S2.Label_; }

      bool operator!=(SiteBasis const& S2) const
      { return Basis_ != S2.Basis_ || Label_ != S2.Label_; }

      void CoerceSymmetryList(SymmetryList const& sl)
      {
	 Basis_.CoerceSymmetryList(sl);
      }

   private:
      typedef std::vector<std::string> LabelType;

      BasisList Basis_;
      pvalue_ptr<std::vector<std::string> > Label_;

   friend PStream::opstream& operator<<(PStream::opstream& out, SiteBasis const& B);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, SiteBasis& B);
   friend void CoerceSymmetryList(SiteBasis& b, SymmetryList const& sl);
   friend SiteBasis adjoint(SiteBasis const& s);
};

std::ostream& operator<<(std::ostream& out, SiteBasis const& Basis);

void show_projections(std::ostream& out, SiteBasis const& Basis);

inline
SiteBasis adjoint(SiteBasis const& s)
{
   SiteBasis Result;
   Result.Basis_ = adjoint(s.Basis_);
   Result.Label_ = s.Label_;
   return Result;
}

//
// Tensor products of site basis
//

class SiteProductBasis
{
   public:
      SiteProductBasis(SiteBasis const& B1, SiteBasis const& B2);

      SiteProductBasis(SiteBasis const& B1, SiteBasis const& B2, 
		       SymmetryList const& ProductSymmetry);

      SiteBasis const& Basis() const { return Basis_; }

      ProductBasis<BasisList, BasisList> const& PBasis() const { return ProductBasis_; }

      std::size_t size() const { return ProductBasis_.size(); }

      SiteBasis Basis1() const { return Basis1_; }
      SiteBasis Basis2() const { return Basis2_; }

   private:
      SiteBasis Basis_, Basis1_, Basis2_;
      ProductBasis<BasisList, BasisList> ProductBasis_;
};

std::ostream& operator<<(std::ostream& out, SiteProductBasis const& Basis);

#endif
