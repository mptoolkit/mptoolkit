// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mpo/product_mpo.h
//
// Copyright (C) 2013-2022 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2024 Jesse Osborne <j.osborne@uqconnect.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER
//
// Matrix product operator for a translationally invariant (or k-dependent)
// product, such as a unitary or non-unitary evolution operator.
//
// Addition isn't defined for ProductMPO.
// We can multiply them and act on an MPS.
//
// Initially they were required to be scalars, but we now allow for a quantum number per unit cell
// (which is the size of the ProductMPO), we can have operators that represent eg products of creation
// operators.

#if !defined(MPTOOLKIT_MPO_PRODUCT_MPO_H)
#define MPTOOLKIT_MPO_PRODUCT_MPO_H

#include "generic_mpo.h"
#include "basic_finite_mpo.h"

class ProductMPO
{
   private:
      typedef GenericMPO data_type;

   public:
      typedef OperatorComponent         value_type;
      typedef data_type::const_iterator const_iterator;
      typedef data_type::iterator       iterator;
      typedef OperatorComponent::basis1_type basis1_type;
      typedef OperatorComponent::basis2_type basis2_type;

      ProductMPO() = default;
      ProductMPO(ProductMPO const& Other) = default;
      ProductMPO(ProductMPO&& Other) = default;

      explicit ProductMPO(int Size) : Data(Size), QShift() {}
      ProductMPO(int Size, QuantumNumbers::QuantumNumber const& Q) : Data(Size), QShift(Q) {}

      // Construction from a generic MPO.  The generic MPO must already be in appropriate form.
      explicit ProductMPO(GenericMPO const& Other)
      : Data(Other), QShift(Other.GetSymmetryList())
      {
         // A ProductMPO constructed in this fashion must be a scalar
         CHECK_EQUAL(Data.Basis1(), Data.Basis2());
      }

      ProductMPO(GenericMPO const& Other, QuantumNumbers::QuantumNumber const& Q)
         : Data(Other), QShift(Q) {}

      ProductMPO& operator=(ProductMPO const& Other) = default;
      ProductMPO& operator=(ProductMPO&& Other) = default;

      // returns the total number of sites this operator contains
      int size() const { return Data.size(); }

      QuantumNumbers::QuantumNumber qshift() const { return QShift; }
      QuantumNumbers::QuantumNumber& qshift() { return QShift; }

      // returns true if this is a zero operator
      bool empty() const { return Data.empty() || Data.front().Basis1().size() == 0; }
      bool is_null() const { return Data.empty() || Data.front().Basis1().size() == 0; }

      // returns true if this MPO is the identity operator, that is, a 1x1 MPO that
      // is a product of identity operators.
      bool is_identity() const;

      // returns true if this is a scalar operator.  In that case, Basis1() == Basis2()
      bool is_scalar() const { return QuantumNumbers::is_scalar(QShift); }

      // returns true if this MPO is a 'string' operator with respect to the unit cell,
      // that is, a 1x1 MPO
      bool is_string() const;

      // returns the left-most basis.
      basis1_type const& Basis1() const { return Data.front().Basis1(); }

      basis2_type const& Basis2() const { return Data.back().Basis2(); }

      value_type& operator[](int n) { return Data[n]; }
      value_type const& operator[](int n) const { return Data[n]; }

      // iterate over the MPOpComponents at each site
      const_iterator begin() const { return Data.begin(); }
      const_iterator end() const { return Data.end(); }

      iterator begin() { return Data.begin(); }
      iterator end() { return Data.end(); }

      value_type& front() { return Data.front(); }
      value_type const& front() const { return Data.front(); }

      value_type& back() { return Data.back(); }
      value_type const& back() const { return Data.back(); }

      // small problem here: if this->is_null(), this will not work.
      SymmetryList GetSymmetryList() const { return Data.front().GetSymmetryList(); }

      // implicit conversion to a const GenericMPO
      operator GenericMPO const&() const { return Data; }

      // returns the list of local hilbert spaces for this operator
      std::vector<BasisList> LocalBasis1List() const
      { return Data.LocalBasis1List(); }
      std::vector<BasisList> LocalBasis2List() const
      { return Data.LocalBasis2List(); }

      // direct access to the GenericMPO
      data_type& data() { return Data; }
      data_type const& data() const { return Data; }

      void check_structure() const;
      void debug_check_structure() const;

      // Make an identity MPO over the given unit cell basis
      static ProductMPO make_identity(std::vector<BasisList> const& Basis);

      // identity operator acting in the given quantum number sector
      static ProductMPO make_identity(std::vector<BasisList> const& Basis,
                                      QuantumNumber const& q);

      static PStream::VersionTag VersionT;

   private:
      data_type Data;
      QuantumNumbers::QuantumNumber QShift;

      friend PStream::opstream& operator<<(PStream::opstream& out, ProductMPO const& op);
      friend PStream::ipstream& operator>>(PStream::ipstream& in, ProductMPO& op);
};

inline
void
ProductMPO::debug_check_structure() const
{
#if !defined(NDEBUG)
   this->check_structure();
#endif
}

// Returns the MPO that is Op1 \otimes Op2.
// PRECONDITION: Op1.Basis2() == Op2.Basis1()
ProductMPO join(ProductMPO const& Op1, ProductMPO const& Op2);

// Repeats Op Count number of times, Op^{\oprod Count}.
// PRECONDITION: Op.Basis2() == Op.Basis1()
ProductMPO repeat(ProductMPO const& Op, int Count);

// multiply by scalar isn't defined for ProductMPO

ProductMPO& operator*=(ProductMPO& x, ProductMPO const& y);

// ProductMPO always represents a scalar, so we don't need the
// 3-argument version of prod()
ProductMPO prod(ProductMPO const& x, ProductMPO const& y);
ProductMPO operator*(ProductMPO const& x, ProductMPO const& y);

inline
ProductMPO dot(ProductMPO const& x, ProductMPO const& y)
{
   return x*y;
}

ProductMPO inner(ProductMPO const& x, ProductMPO const& y);

ProductMPO outer(ProductMPO const& x, ProductMPO const& y);

// power of an operator.  Requires n > 1.
ProductMPO pow(ProductMPO const& x, int n);

// Conjugate
ProductMPO conj(ProductMPO const& x);

// Adjoint
ProductMPO adjoint(ProductMPO const& x);
ProductMPO inv_adjoint(ProductMPO const& x);

// constructs a ProductMPO as a repeated string of some BasicFiniteMPO
ProductMPO string(BasicFiniteMPO const& Op);

ProductMPO gauge_flip(ProductMPO const& Op);

// Constructs a ProductMPO as the infinite product of translations of Op.
// Op.size() must be an integer multiple of UnitCellSize,
// and the local basis of Op[i] must be the same as the local basis of Op[i%UnitCellSize]
// This does the products in the left-to-right sense, that is, if our operator is
// AB then the product is ... * A(0)B(1) * A(1)B(2) * A(2)B(3) * ...
// = B(0)A(0) B(1)A(1) B(2)A(2) ...
ProductMPO prod_unit_left_to_right(BasicFiniteMPO const& Op, int UnitCellSize);


inline
ProductMPO prod_unit(BasicFiniteMPO const& Op, int UnitCellSize)
{
   return prod_unit_left_to_right(Op, UnitCellSize);
}

// Constructs a ProductMPO as the infinite product of translations of Op.
// Op.size() must be an integer multiple of UnitCellSize,
// and the local basis of Op[i] must be the same as the local basis of Op[i%UnitCellSize]
// This does the products in the right-to-left sense, that is, if our operator is
// AB then the product is .... * A(2)B(3) * A(1)B(2) * A(0)B(1) * ...
// = A(0)B(0) A(1)B(1) A(2)B(2) ...
ProductMPO prod_unit_right_to_left(BasicFiniteMPO const& Op, int UnitCellSize);

// Construct the operator to do a right translation by one site, given the
// array of local basis states
ProductMPO translate_right(std::vector<BasisList> const& LocalBasis);

// optimize the representation - no idea how to do this!
//void optimize(ProductMPO& Op);

// Does a N-1 coarse graining of an operator.  The length must be a multiple of N
inline
ProductMPO coarse_grain(ProductMPO const& Op, int N)
{
   return ProductMPO(coarse_grain(Op.data(), N), Op.qshift());
}

// output to a stream
std::ostream& operator<<(std::ostream& out, ProductMPO const& x);

// prints the structure of the component, as an 'x' for a non-zero
// component or blank for a zero component
void print_structure(ProductMPO const& Op, std::ostream& out, double UnityEpsilon);

inline
void print_structure(ProductMPO const& Op, std::ostream& out)
{
   print_structure(Op, out, DefaultClassifyUnityEpsilon);
}

#endif
