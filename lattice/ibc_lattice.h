// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// lattice/ibc_lattice.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
// -*- C++ -*- $Id$

// An InfiniteLattice is defined over a UnitCell, and also defines operators
// with infinite support (ie Triangular mpo's).
//
// We really should also define also a ProductMPO.  This would represent the
// infinite product of operators defined on the unit cell, and also have a QShift, with
// the requirement that Basis1() == delta_shift(Basis2(), QShift).
//
// The unit cell of the operators is allowed to a a multiple of the lattice UnitCell.

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

      InfiniteLattice();

      explicit InfiniteLattice(UnitCell const& uc);

      InfiniteLattice(std::string const& Description, UnitCell const& uc);

      UnitCell const& GetUnitCell() const { return UnitCell_; }

      std::string description() const { return Description_; }
      void set_description(std::string s) { Description_ = s; }

      // Set the command_line used to construct the lattice.  Also sets the timestamp.
      void set_command_line(int argc, char** argv);

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

      // returns true if the given operator exists, and is a TriangularMPO
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
      TriangularMPO const& as_triangular_mpo(std::string const& Op) const;
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
      std::string CommandLine_;
      std::string Timestamp_;
      UnitCell UnitCell_;
      OperatorListType Operators_;
      ArgumentListType Arguments_;
      FunctionListType Functions_;

   friend PStream::opstream& operator<<(PStream::opstream& out, InfiniteLattice const& L);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, InfiniteLattice& L);
};

// Constucts a TriangularMPO from the summation over unit cell translations of a finite MPO.
// The Op must have a size() that is a multiple of SiteListTypeSize, which must itself be an
// integer multiple of SiteListType.size().
// The aux basis for JW is assumed to be compatible with Op -- that is, JW.qn2() == Op.qn1()
TriangularMPO sum_unit(SiteListType const& SiteList, FiniteMPO const& JW, FiniteMPO const& Op, int UnitCellSize);

TriangularMPO sum_unit(SiteListType const& SiteList, FiniteMPO const& Op, LatticeCommute Com, int UnitCellSize);

TriangularMPO sum_unit(SiteListType const& SiteList, FiniteMPO const& Op, LatticeCommute Com);

TriangularMPO sum_unit(UnitCellMPO const& Op, int UnitCellSize);

TriangularMPO sum_unit(UnitCellMPO const& Op);

// Variant of sum_unit where we add the kink operator (generally will be unitary) to the left hand side
TriangularMPO sum_kink(SiteListType const& SiteList, FiniteMPO const& Kink,
		       FiniteMPO const& Op, LatticeCommute Com, int UnitCellSize);

TriangularMPO sum_kink(UnitCellMPO const& Kink, UnitCellMPO const& Op, int UnitCellSize);

TriangularMPO sum_kink(UnitCellMPO const& Kink, UnitCellMPO const& Op);

// sum_k

TriangularMPO sum_k(SiteListType const& SiteList, std::complex<double> const& k,
		       FiniteMPO const& Op, LatticeCommute Com, int UnitCellSize);

TriangularMPO sum_k(std::complex<double> const& k, UnitCellMPO const& Op, int UnitCellSize);

TriangularMPO sum_k(std::complex<double> const& k, UnitCellMPO const& Op);

// sum_string
// Constsructs a TriangularMPO of the form
// A(0)*C(1) + A(0)*B(1)*C(2) + A(0)*B(1)*B(2)*C(3) + ....
// plus translations.
// PRECONDITION(JW.qn2() == Op1.qn1() && Op1.qn2() == String.qn1() && String.qn2() == Op1.qn1()
// && is_scalar(Op2.qn2()))
TriangularMPO sum_string(SiteListType const& SiteList, FiniteMPO const& JW, FiniteMPO const& Op1,
			 FiniteMPO const& String, FiniteMPO const& Op2, int UnitCellSize,
			 QuantumNumbers::QuantumNumber q);

// This version of sum_string takes UnitCellMPO's for the operator arguments.  The String term
// must be a scalar with bosonic commutation, and cannot be any longer than UnitCellSize.
TriangularMPO sum_string(UnitCellMPO const& Op1_, UnitCellMPO const& String_, UnitCellMPO const& Op2_, 
			 int UnitCellSize,
			 QuantumNumbers::QuantumNumber q);

TriangularMPO sum_string_dot(UnitCellMPO const& Op1_, UnitCellMPO const& String_, UnitCellMPO const& Op2_, 
			     int UnitCellSize);

TriangularMPO sum_string_dot(UnitCellMPO const& Op1_, UnitCellMPO const& String_, UnitCellMPO const& Op2_);

TriangularMPO sum_string_inner(UnitCellMPO const& Op1_, UnitCellMPO const& String_, UnitCellMPO const& Op2_, 
			       int UnitCellSize);

TriangularMPO sum_string_inner(UnitCellMPO const& Op1_, UnitCellMPO const& String_, UnitCellMPO const& Op2_);

// Constructs a zero triangular MPO
TriangularMPO make_zero(SiteListType const& SiteList);

#endif
