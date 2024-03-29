// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// lattice/infinitelattice.h
//
// Copyright (C) 2014-2022 Ian McCulloch <ian@qusim.net>
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

// An InfiniteLattice is defined over a UnitCell, and also defines operators
// with infinite support (ie Triangular mpo's).
//
// We really should also define also a ProductMPO.  This would represent the
// infinite product of operators defined on the unit cell, and also have a QShift, with
// the requirement that Basis1() == delta_shift(Basis2(), QShift).
//
// The unit cell of the operators is allowed to a multiple of the lattice UnitCell.
// We do some fancy tracking to keep the UnitCell by pointer where possible

#if !defined(MPTOOLKIT_LATTICE_INFINITELATTICE_H)
#define MPTOOLKIT_LATTICE_INFINITELATTICE_H

#include "unitcell.h"
#include "mpo/infinite_mpo.h"
#include <vector>
#include "pheap/pvalueptr.h"
#include "lattice/function.h"
#include "lattice/operator_descriptions.h"

// Lattice version number for streaming
extern PStream::VersionTag LatticeVersion;

class InfiniteLattice
{
   public:
      typedef InfiniteMPO                operator_type;
      typedef Function::OperatorFunction function_type;
      typedef Function::Argument         argument_type;

   private:
      typedef std::map<std::string, operator_type>    OperatorListType;
      typedef std::map<std::string, function_type>    FunctionListType;
      typedef Function::ArgumentList                  ArgumentListType;

   public:
      typedef OperatorListType::iterator        operator_iterator;
      typedef OperatorListType::const_iterator  const_operator_iterator;

      typedef FunctionListType::iterator        function_iterator;
      typedef FunctionListType::const_iterator  const_function_iterator;

      typedef ArgumentListType::iterator        argument_iterator;
      typedef ArgumentListType::const_iterator  const_argument_iterator;

      typedef std::vector<std::pair<std::string, std::string>> authors_type;

      InfiniteLattice();

      explicit InfiniteLattice(UnitCell* uc);
      explicit InfiniteLattice(UnitCell&& uc);
   //      explicit InfiniteLattice(UnitCell const& uc);  // maybe confusing to have this?

      InfiniteLattice(std::string const& Description, UnitCell& uc);
      InfiniteLattice(std::string const& Description, UnitCell&& uc);

      InfiniteLattice(InfiniteLattice const& Other);
      InfiniteLattice(InfiniteLattice&& Other);

      InfiniteLattice& operator=(InfiniteLattice const& Other);
      InfiniteLattice& operator=(InfiniteLattice&& Other);

      ~InfiniteLattice();

      UnitCell& GetUnitCell() { return *UnitCell_; }
      UnitCell const& GetUnitCell() const { return *UnitCell_; }

      std::string description() const { return Description_; }
      void set_description(std::string s) { Description_ = s; }

      // Set the command_line used to construct the lattice.  Also sets the timestamp.
      void set_command_line(int argc, char** argv);

      // lattice file author information
      authors_type& authors() { return Authors_; }
      authors_type const& authors() const { return Authors_; }

      std::string command_line() const { return CommandLine_; }
      std::string timestamp() const { return Timestamp_; }

      SymmetryList GetSymmetryList() const { return this->GetUnitCell().GetSymmetryList(); }

      // operators

      const_operator_iterator begin_operator() const
      { return Operators_.begin(); }
      const_operator_iterator end_operator() const
      { return Operators_.end(); }

      // returns true if the given operator exits
      bool operator_exists(std::string const& s) const;

      // returns true if the given operator exists, and is a BasicTriangularMPO
      bool triangular_operator_exists(std::string const& s) const;

      // returns true if the given operator exists, and is a ProductMPO
      bool product_operator_exists(std::string const& s) const;

      // Sets the description of the operators according to Desc.  Prints a warning to cerr
      // if there are extra operators named in Desc that are not part of the lattice,
      // or if any operators don't have a name.
      void set_operator_descriptions(OperatorDescriptions const& Desc);

      // Lookup the named operator
      InfiniteMPO& operator[](std::string const& Op);
      InfiniteMPO const& operator[](std::string const& Op) const;
      BasicTriangularMPO const& as_basic_triangular_mpo(std::string const& Op) const;
      ProductMPO const& as_product_mpo(std::string const& Op) const;

      // arguments

      bool arg_empty() const { return Arguments_.empty(); }

      const_argument_iterator begin_arg() const { return Arguments_.begin(); }
      const_argument_iterator end_arg() const { return Arguments_.end(); }

      const_argument_iterator find_arg(std::string const& s) const
      { return Arguments_.find(s); }

      argument_iterator find_arg(std::string const& s)
      { return Arguments_.find(s); }

      std::complex<double> arg(std::string const& a) const;

      std::complex<double>& arg(std::string const& a)
      { return Arguments_[a]; }

      ArgumentListType const& args() const { return Arguments_; }

      // functions

      bool function_empty() const { return Functions_.empty(); }

      // iterators
      const_function_iterator begin_function() const { return Functions_.begin(); }
      const_function_iterator end_function() const { return Functions_.end(); }

      // find
      const_function_iterator find_function(std::string const& s) const
      { return Functions_.find(s); }

      bool function_exists(std::string const& s) const
      { return Functions_.find(s) != Functions_.end(); }

      // lookup
      function_type& func(std::string const& s) { return Functions_[s]; }

      function_type func(std::string const& s) const;

      // evaluate a function
      InfiniteMPO
      eval_function(Function::OperatorFunction const& Func,
                    Function::ParameterList const& Params) const;

      InfiniteMPO
      eval_function(std::string const& Func,
                    Function::ParameterList const& Params) const;

   private:
      std::string Description_;
      authors_type Authors_;
      std::string CommandLine_;
      std::string Timestamp_;
      UnitCell* UnitCell_;
      bool OwnUnitCell_;               // true if UnitCell_ is owned by us on the heap
      OperatorListType Operators_;
      ArgumentListType Arguments_;
      FunctionListType Functions_;

   friend PStream::opstream& operator<<(PStream::opstream& out, InfiniteLattice const& L);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, InfiniteLattice& L);
};

// Constructs a zero triangular MPO
BasicTriangularMPO make_zero(SiteListType const& SiteList);

// Constucts a BasicTriangularMPO from the summation over unit cell translations of a finite MPO.
// This is a special case of sum_kink() where the kink operator is the identity.
// The Op must have a size() that is a multiple of SiteListTypeSize, which must itself be an
// integer multiple of SiteListType.size().
// The aux basis for JW is assumed to be compatible with Op -- that is, JW.qn2() == Op.qn1()
BasicTriangularMPO sum_unit(BasicFiniteMPO const& JW, BasicFiniteMPO const& Op, int UnitCellSize);

BasicTriangularMPO sum_unit(UnitCellMPO const& Op, int UnitCellSize);

BasicTriangularMPO sum_unit(UnitCellMPO const& Op);

BasicTriangularMPO make_triangular(UnitCellMPO const& Op);

inline
ProductMPO prod_unit_left_to_right(UnitCellMPO const& Op, int UnitCellSize)
{
   CHECK(UnitCellSize % Op.coarse_grain_factor() == 0);
   return prod_unit_left_to_right(Op.MPO(), UnitCellSize / Op.coarse_grain_factor());
}

inline
ProductMPO prod_unit_right_to_left(UnitCellMPO const& Op, int UnitCellSize)
{
   CHECK(UnitCellSize % Op.coarse_grain_factor() == 0);
   return prod_unit_right_to_left(Op.MPO(), UnitCellSize / Op.coarse_grain_factor());
}

inline
ProductMPO prod_unit(UnitCellMPO const& Op)
{
   return prod_unit_left_to_right(Op, Op.GetSiteList()->size());
}

inline
ProductMPO prod_unit(UnitCellMPO const& Op, int UnitCellSize)
{
   return prod_unit_left_to_right(Op, UnitCellSize);
}

// Variant of sum_unit where we add the kink operator (generally will be unitary) to the left hand side
BasicTriangularMPO sum_kink(UnitCellMPO const& Kink, UnitCellMPO const& Op, int UnitCellSize);

BasicTriangularMPO sum_kink(UnitCellMPO const& Kink, UnitCellMPO const& Op);

// sum_k
// Consruct an operator at finite momentum.  This is a special case of sum_kink where the kink operator is
// e^{ik} times the identity.
BasicTriangularMPO sum_k(std::complex<double> const& k, UnitCellMPO const& Op, int UnitCellSize);

BasicTriangularMPO sum_k(std::complex<double> const& k, UnitCellMPO const& Op);

// sum_string
// Constsructs a BasicTriangularMPO of the form
// A*C + A*B*C + A*B*B*C + ...
// plus translations.
// PRECONDITION(JW.qn2() == Op1.qn1() && Op1.qn2() == String.qn1() && String.qn2() == Op1.qn1()
// && is_scalar(Op2.qn2()))
BasicTriangularMPO sum_string(SiteListType const& SiteList, BasicFiniteMPO const& JW, BasicFiniteMPO const& Op1,
                         BasicFiniteMPO const& String, BasicFiniteMPO const& Op2, int UnitCellSize,
                         QuantumNumbers::QuantumNumber q);

// This version of sum_string takes UnitCellMPO's for the operator arguments.  The String term
// must be a scalar with bosonic commutation, and cannot be any longer than UnitCellSize.
BasicTriangularMPO sum_string(UnitCellMPO const& Op1_, UnitCellMPO const& String_, UnitCellMPO const& Op2_,
                         int UnitCellSize,
                         QuantumNumbers::QuantumNumber q);

BasicTriangularMPO sum_string(UnitCellMPO const& Op1_, UnitCellMPO const& String_, UnitCellMPO const& Op2_);

// 5-term variation of sum_string
BasicTriangularMPO sum_string(SiteListType const& SiteList, BasicFiniteMPO const& JW, BasicFiniteMPO const& Op1,
                         BasicFiniteMPO const& String1, BasicFiniteMPO const& Op2,
                         BasicFiniteMPO const& String2, BasicFiniteMPO const& Op3,
                         int UnitCellSize, QuantumNumbers::QuantumNumber q);

BasicTriangularMPO sum_string(UnitCellMPO const& Op1_, UnitCellMPO const& String1_, UnitCellMPO const& Op2_,
                              UnitCellMPO const& String2_, UnitCellMPO const& Op3_,
                              int UnitCellSize,
                              QuantumNumbers::QuantumNumber q);

BasicTriangularMPO sum_string(UnitCellMPO const& Op1_, UnitCellMPO const& String1_, UnitCellMPO const& Op2_,
                              UnitCellMPO const& String2_, UnitCellMPO const& Op3_);

// sum_string_dot
// variant of sum_string where we take the (non-abelian) dot product between operators A and C.
BasicTriangularMPO sum_string_dot(UnitCellMPO const& Op1_, UnitCellMPO const& String_, UnitCellMPO const& Op2_,
                             int UnitCellSize);

BasicTriangularMPO sum_string_dot(UnitCellMPO const& Op1_, UnitCellMPO const& String_, UnitCellMPO const& Op2_);

// sum_string_inner
// variant of sum_string where we take the (non-abelian) inner product between operators A and C.
BasicTriangularMPO sum_string_inner(UnitCellMPO const& Op1_, UnitCellMPO const& String_, UnitCellMPO const& Op2_,
                               int UnitCellSize);

BasicTriangularMPO sum_string_inner(UnitCellMPO const& Op1_, UnitCellMPO const& String_, UnitCellMPO const& Op2_);

// sum_partial
// Constructs a new MPO that is all translations of partial sums of an existing triangular MPO.  Eg,
// let M = sum_i M(i) be a triangular MPO, where M(j) means all terms that have support on sites up to j (inclusive).
// Then sum_partial(M, P) is
// sum_j (sum_{i < j} M(i)) P(j)
//
// Concrete example: sum_partial(sum_unit(Sx), Sy) is the two-dimensional sum_{i<j} Sx(i) Sy(j)
// which is equivalent to sum_string(Sx, I, Sy)
BasicTriangularMPO sum_partial(BasicTriangularMPO const& Op, UnitCellMPO const& Pivot, int UnitCellSize);

BasicTriangularMPO sum_partial(BasicTriangularMPO const& Op, UnitCellMPO const& Pivot);

#endif
