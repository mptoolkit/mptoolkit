// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// lattice/latticesite.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER
/* -*- C++ -*- $Id$

  Created 2004-09-05 Ian McCulloch

  originally a 'site block' in DMRG terminology, now a LatticeSite is
  a collection of SiteOperator's defined on some fixed Hilbert space.

  This also has a description, which is just used for informational
  purposes, eg "Fermion site", "Boson site", etc
*/

#if !defined(MPTOOLKIT_LATTICE_LATTICESITE_H)
#define MPTOOLKIT_LATTICE_LATTICESITE_H

#include "pstream/pstream.h"
#include "pheap/pvalueptr.h"
#include "quantumnumbers/symmetrylist.h"
#include "siteoperator.h"
#include "function.h"
#include "operator_descriptions.h"
#include <boost/variant.hpp>
#include <map>

using QuantumNumbers::SymmetryList;

class LatticeSite
{
   public:
      typedef SiteOperator                operator_type;
      typedef Function::OperatorFunction  function_type;
      typedef Function::Argument         argument_type;

   private:
      typedef std::map<std::string, operator_type>    OperatorListType;
      typedef std::map<std::string, function_type>    FunctionListType;
      typedef Function::ArgumentList                  ArgumentListType;

   public:
      typedef operator_type::basis1_type        basis1_type;
      typedef operator_type::basis2_type        basis2_type;
      typedef OperatorListType::iterator        operator_iterator;
      typedef OperatorListType::const_iterator  const_operator_iterator;

      typedef FunctionListType::iterator        function_iterator;
      typedef FunctionListType::const_iterator  const_function_iterator;

      typedef ArgumentListType::iterator        argument_iterator;
      typedef ArgumentListType::const_iterator  const_argument_iterator;

      LatticeSite() : pImpl(new ImplType()) {}

      explicit LatticeSite(std::string const& Description) 
	 : pImpl(new ImplType(Description)) {}

      LatticeSite(std::string const& Description, ArgumentListType const& Args) 
	 : pImpl(new ImplType(Description, Args)) {}

      std::string const& Description() const { return pImpl->Description; }
      void SetDescription(std::string const& s) { pImpl.mutate()->Description = s; }

      // precondition: !empty()
      SymmetryList GetSymmetryList() const;

      // precondition: !empty()
      basis1_type const& Basis1() const;
      basis2_type const& Basis2() const;

      // Sets the description of the operators according to Desc.  Prints a warning to cerr
      // if there are extra operators named in Desc that are not part of the lattice,
      // or if any operators don't have a name.
      void set_operator_descriptions(OperatorDescriptions const& Desc);

      // operators
      // returns the identity operator
      operator_type identity() const;

      // operator_empty() function is redundant, since we always have at least the identity operator

      const_operator_iterator begin_operator() const { return pImpl->Operators.begin(); }
      const_operator_iterator end_operator() const { return pImpl->Operators.end(); }

      const_operator_iterator find_operator(std::string const& s) const 
      { return pImpl->Operators.find(s); }

      bool operator_exists(std::string const& s) const 
      { return pImpl->Operators.find(s) != pImpl->Operators.end(); }

      operator_type& operator[](std::string const& s) { return pImpl.mutate()->Operators[s]; }
      operator_type const& operator[](std::string const& s) const;

      // arguments

      bool arg_empty() const { return pImpl->Arguments.empty(); }
      
      const_argument_iterator begin_arg() const { return pImpl->Arguments.begin(); }
      const_argument_iterator end_arg() const { return pImpl->Arguments.end(); }

      const_argument_iterator find_arg(std::string const& s) const
      { return pImpl->Arguments.find(s); }

      argument_iterator find_arg(std::string const& s)
      { return pImpl.mutate()->Arguments.find(s); }

      std::complex<double> arg(std::string const& a) const;

      std::complex<double>& arg(std::string const& a) 
      { return pImpl.mutate()->Arguments[a]; }

      // functions

      bool function_empty() const { return pImpl->Functions.empty(); }

      const_function_iterator begin_function() const { return pImpl->Functions.begin(); }
      const_function_iterator end_function() const { return pImpl->Functions.end(); }

      const_function_iterator find_function(std::string const& s) const 
      { return pImpl->Functions.find(s); }

      bool function_exists(std::string const& s) const 
      { return pImpl->Functions.find(s) != pImpl->Functions.end(); }

      function_type& func(std::string const& s) { return pImpl.mutate()->Functions[s]; }

      function_type const& func(std::string const& s) const;

      // evaluate a function
      boost::variant<operator_type, std::complex<double> >
      eval_function(Function::OperatorFunction const& Func, 
		    Function::ParameterList const& Params) const;

      boost::variant<operator_type, std::complex<double> >
      eval_function(std::string const& Func, 
		    Function::ParameterList const& Params) const;

      void check_structure() const;

      void debug_check_structure() const;

      struct ImplType
      {
         std::string Description;
         OperatorListType Operators;
	 ArgumentListType Arguments;
	 FunctionListType Functions;

	 ImplType() {}
	 ImplType(std::string const& Desc_) : Description(Desc_) {}
	 ImplType(std::string const& Desc_,
		  Function::ArgumentList const& Args)
	    : Description(Desc_), Arguments(Args) {}

         friend PStream::opstream& operator<<(PStream::opstream& out, ImplType const& Impl);
         friend PStream::ipstream& operator>>(PStream::ipstream& in, ImplType& Impl);
      };

   private:
      typedef pvalue_ptr<ImplType> ptr_type;
      ptr_type pImpl;

   friend PStream::opstream& operator<<(PStream::opstream& out, LatticeSite const& B);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, LatticeSite& B);

   friend void CoerceSymmetryListInPlace(LatticeSite& s, QuantumNumbers::SymmetryList const& sl);
   friend LatticeSite CoerceSymmetryList(LatticeSite const& s, QuantumNumbers::SymmetryList const& sl);
};

std::ostream& operator<<(std::ostream& out, LatticeSite const& s);

using Function::par;  // shortcut for constructing parameters

// proxy function for an operator
class SiteOperatorFunction
{
   public:
      typedef std::complex<double> complex;

      SiteOperatorFunction(LatticeSite const* Site, std::string const& Name,
			   Function::OperatorFunction const& Func);

      SiteOperator operator()(Function::Parameter const& arg1);
      SiteOperator operator()(Function::Parameter const& arg1, 
			      Function::Parameter const& arg2);
      SiteOperator operator()(Function::Parameter const& arg1, 
			      Function::Parameter const& arg2,
			      Function::Parameter const& arg3);
      SiteOperator operator()(Function::Parameter const& arg1, 
			      Function::Parameter const& arg2,
			      Function::Parameter const& arg3,
			      Function::Parameter const& arg4);
      SiteOperator operator()(Function::Parameter const& arg1, 
			      Function::Parameter const& arg2,
			      Function::Parameter const& arg3,
			      Function::Parameter const& arg4,
			      Function::Parameter const& arg5);
      SiteOperator operator()(Function::Parameter const& arg1, 
			      Function::Parameter const& arg2,
			      Function::Parameter const& arg3,
			      Function::Parameter const& arg4,
			      Function::Parameter const& arg5,
			      Function::Parameter const& arg6);
};

// This is used by UnitCell and UnitCellMPO
typedef std::vector<LatticeSite> SiteListType;
typedef pvalue_ptr<SiteListType> SiteListPtrType;

// utility to gt a vector of the local basis from a SiteListType
std::vector<BasisList>
Basis1FromSiteList(SiteListType const& s);

#endif
