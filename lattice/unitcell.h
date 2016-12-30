// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// lattice/unitcell.h
//
// Copyright (C) 2014-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

// A UnitCell is an array of LatticeSite's, that represents a primitive unit for a lattice.
// The sites in the lattice are denoted [0], [1], .... using square brackets, and zero-based counting.
//
// For translationally invariant systems, this is the basic unit that defines the repeated unit.
//
// the unit cell defines operators with finite support.  Eg, the UnitCell might contain N LatticeSite's,
// each of which has a spin operator S.  We refer to the n'th LatticeSite operator as S[n].
// Typically a UnitCell will also define its own operators that act on the UnitCell as a whole,
// for example S = S[0]+S[1]+...+S[N].
// Once we get to a Lattice, as a collection of UnitCell's, we would refer to the total spin operator at
// the i'th UnitCell as S(i).  An individual spin at the n'th site of the i'th UnitCell would be referred
// to as S(i)[n].
//
// The unit cell of the operators is allowed to be a multiple of the lattice UnitCell. In this case,
// we would refer to a specific operator on the lattice as the left-most unit cell where
// it has non-trivial support.  This lets us define UnitCell operators that span more than one unit cell,
// eg for bond operators.
//
// The UnitCell always contains the global identity operator "I".  This is added by the UnitCell
// constructor.
//
// If the UnitCell is a single site, we copy all of the operators, arguments and functions
// to act on the UnitCell also.  (This replaces the older, flawed, method of delegating to
// the SiteOperator in the special case of size()==1, which is broken because, eg, iterators
// don't work properly.)
// **** There is a problem here, when we lift site operators to a unit cell, if we parse
// a function of the form U(site){parameters}, then we don't have the site index when parsing.
// we need to parse the operator in the context of a SiteOperator.
//
// If UnitCell's are join()'ed to make a new unit cell, then operators, arguments, and functions
// are NOT carried over to the joined UnitCell.  There is no unique way of defining what such
// objects should look like.

#if !defined(MPTOOLKIT_LATTICE_UNITCELL_H)
#define MPTOOLKIT_LATTICE_UNITCELL_H

#include "lattice/latticesite.h"
#include "lattice/unitcell_mpo.h"
#include <vector>
#include <map>
#include <string>

class UnitCell
{
   public:
      typedef LatticeSite               value_type;
      typedef std::vector<LatticeSite>  data_type;
      typedef data_type::const_iterator const_iterator;

      typedef UnitCellMPO                operator_type;
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

      UnitCell();

      // Copy ctor is non-trivial since we need to reset the back-pointer in the UnitCellMPO's
      UnitCell(UnitCell const& Other);

      UnitCell(LatticeSite const& s);
      UnitCell(LatticeSite const& s, LatticeSite const& t);
      UnitCell(LatticeSite const& s, LatticeSite const& t, LatticeSite const& u);
      UnitCell(LatticeSite const& s, LatticeSite const& t, LatticeSite const& u, LatticeSite const& v);

      UnitCell(std::vector<LatticeSite> const& s);

      UnitCell(SymmetryList const& sl, LatticeSite const& s);
      UnitCell(SymmetryList const& sl, LatticeSite const& s, LatticeSite const& t);
      UnitCell(SymmetryList const& sl, LatticeSite const& s, LatticeSite const& t,
               LatticeSite const& u);
      UnitCell(SymmetryList const& sl, LatticeSite const& s, LatticeSite const& t,
               LatticeSite const& u, LatticeSite const& v);

      // Construct by joining unit cells.
      // NOTE: operators, arguments and functions are NOT inherited when
      // joining UnitCell's

      UnitCell(int RepeatCount, UnitCell const& l);
      UnitCell(UnitCell const& x1, UnitCell const& x2);
      UnitCell(UnitCell const& x1, UnitCell const& x2, UnitCell const& x3);
      UnitCell(UnitCell const& x1, UnitCell const& x2, UnitCell const& x3, UnitCell const& x4);

      UnitCell(int Size, LatticeSite const& s);

      UnitCell& operator=(UnitCell const& Other);

      SymmetryList GetSymmetryList() const { return Sites->front().GetSymmetryList(); }

      // LatticeSite methods

      bool empty() const { return Sites->empty(); }
      int size() const { return Sites->size(); }
      value_type const& front() const { return Sites->front(); }
      value_type const& back() const { return Sites->back(); }

      const_iterator begin() const { return Sites->begin(); }
      const_iterator end() const { return Sites->end(); }

      // returns the i'th LatticeSite
      LatticeSite const& operator[](int i) const;

#if 0
      // this turns out to be not needed so far.
      // Visitor pattern tools for iterating over the LatticeSite's
      template <typename Visitor>
      typename Visitor::result_type
      apply_for_each_site(Visitor const& v) const;

      template <typename Visitor>
      typename Visitor::result_type
      apply_for_each_site(Visitor const& v);
#endif

      SiteListPtrType const& GetSiteList() const { return Sites; }

      // returns the local basis at site i of the lattice
      SiteBasis LocalBasis(int i) const
      { return this->operator[](i).Basis2(); }

      // operators

      // iterators
      const_operator_iterator begin_operator() const { return Operators.begin(); }
      const_operator_iterator end_operator() const { return Operators.end(); }

      // returns the identity operator on the UnitCell
      operator_type identity() const;

      // returns true if the named operator exists
      bool operator_exists(std::string const& s) const;

      // lookup a unit cell operator
      operator_type as_operator(std::string const& Op) const;

      // lookup a unit cell operator,
      // adds it if it doesn't already exist.
      operator_type& operator[](std::string const& Op);

      operator_type operator[](std::string const& Op) const;

      // Assigns the operator, given an offset.  Eg, Name(0) = Op
      void assign_operator(std::string const& Name, operator_type Op, int Offset = 0);

      // returns the operator, shifted to to the given cell number
      operator_type operator()(std::string const& Op, int Cell) const;

      // returns true if the specified operator Op[n] exists at the given Site of the UnitCell
      bool local_operator_exists(std::string const& Op, int Site) const;

      // lookup a local operator, at cell 0
      operator_type local_operator(std::string const& Op, int Site) const;

      // Lookup a local operator on a given unit cell index
      operator_type local_operator(std::string const& Op, int Cell, int Site) const;

      // Given a SiteOperator, map it into a UnitCellMPO at the given (Cell)[Site] coordinates
      operator_type map_local_operator(SiteOperator const& Op, int Cell, int Site) const;

      // arguments

      bool arg_empty() const { return Arguments.empty(); }

      const_argument_iterator begin_arg() const { return Arguments.begin(); }
      const_argument_iterator end_arg() const { return Arguments.end(); }

      const_argument_iterator find_arg(std::string const& s) const
      { return Arguments.find(s); }

      argument_iterator find_arg(std::string const& s)
      { return Arguments.find(s); }

      std::complex<double> arg(std::string const& a) const;

      std::complex<double>& arg(std::string const& a)
      { return Arguments[a]; }

      // functions

      bool function_empty() const { return Functions.empty(); }

      // iterators
      const_function_iterator begin_function() const { return Functions.begin(); }
      const_function_iterator end_function() const { return Functions.end(); }

      // find
      const_function_iterator find_function(std::string const& s) const
      { return Functions.find(s); }

      bool function_exists(std::string const& s) const
      { return Functions.find(s) != Functions.end(); }

      // lookup
      function_type& func(std::string const& s) { return Functions[s]; }

      function_type func(std::string const& s) const;

      // evaluate a function
      boost::variant<operator_type, std::complex<double> >
      eval_function(Function::OperatorFunction const& Func, int Cell,
                    Function::ParameterList const& Params) const;

      // Evaluate a function.  If Func is not a UnitCell function and
      // the UnitCell.size() == 1, then evaluate Func as a SiteOperator function.
      boost::variant<operator_type, std::complex<double> >
      eval_function(std::string const& Func, int Cell,
                    Function::ParameterList const& Params) const;

      boost::variant<operator_type, std::complex<double> >
      eval_local_function(std::string const& Func, int Cell, int Site,
                          Function::ParameterList const& Params) const;

      // special UnitCell operators

      // Returns an MPO that effects a swap gate between sites i and j
      operator_type swap_gate(int i, int j) const;

      // Returns an MPO that effects a swap gate between different unit cells
      operator_type swap_gate(int iCell, int i, int jCell, int j) const;

      // Returns an MPO that effects a swap gate between sites i and j,
      // this version ignores fermions / JW strings
      operator_type swap_gate_no_sign(int i, int j) const;

      // Returns an MPO that effects a swap gate between different unit cells
      // this version ignores fermions / JW strings
      operator_type swap_gate_no_sign(int iCell, int i, int jCell, int j) const;

      // returns the commutation attribute of the operator, equivalent
      // to as_operator(OpName, n).Commute()
      LatticeCommute Commute(std::string const& OpName, int n) const;

      void debug_check_structure() const;
      void check_structure() const;

   private:
      // Set the default operators - currently the Identity operator I.
      // Called after Sites is initialized.
      void SetDefaultOperators();

      pvalue_ptr<SiteListType> Sites;
      OperatorListType Operators;
      ArgumentListType Arguments;
      FunctionListType Functions;


      friend PStream::opstream& operator<<(PStream::opstream& out, UnitCell const& L);
      friend PStream::ipstream& operator>>(PStream::ipstream& in, UnitCell& L);
};

// Make a new lattice out of RepeatCount copies of x
// NOTE: operators, arguments and functions are NOT inherited when
// joining UnitCell's
UnitCell repeat(UnitCell const& x, int RepeatCount);

// Make a new lattice out of the join of x and y
// NOTE: operators, arguments and functions are NOT inherited when
// joining UnitCell's
UnitCell join(UnitCell const& x, UnitCell const& y);
UnitCell join(UnitCell const& x, UnitCell const& y, UnitCell const& z);
UnitCell join(UnitCell const& x, UnitCell const& y, UnitCell const& z, UnitCell const& w);
UnitCell join(UnitCell const& x, UnitCell const& y, UnitCell const& z, UnitCell const& w,
              UnitCell const& v);

#if 0 // These operators probably don't make much sense
bool
operator==(UnitCell const& u1, UnitCell const& u2);

bool
operator!=(UnitCell const& u1, UnitCell const& u2);
#endif

std::ostream&
operator<<(std::ostream& out, UnitCell const& u);

#include "unitcell.cc"

#endif
