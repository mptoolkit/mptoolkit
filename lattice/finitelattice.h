// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// lattice/finitelattice.h
//
// Copyright (C) 2014-2017 Ian McCulloch <ian@qusim.net>
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

// A FiniteLattice is an array of Unit Cells.  This uses run-length compression to handle
// large arrays.

#if !defined(MPTOOLKIT_LATTICE_FINITELATTICE_H)
#define MPTOOLKIT_LATTICE_FINITELATTICE_H

#include <vector>
#include "quantumnumbers/braidgroup.h"
#include "finitelatticempo.h"
#include "unitcell.h"
#include "mpo/basic_finite_mpo.h"
#include "pheap/pvalueptr.h"
#include "lattice/function.h"
#include "lattice/operator_descriptions.h"

// Lattice version number for streaming
extern PStream::VersionTag FiniteLatticeVersion;

typedef run_length_compressed<UnitCell> UnitCellListType;
typedef pvalue_ptr<UnitCellListType> UnitCellPtrType;

class FiniteLattice
{
   public:
      typedef FiniteLatticeMPO           operator_type;
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

      FiniteLattice();

      explicit FiniteLattice(UnitCell const& uc);

      FiniteLattice(std::string const& Description);
      FiniteLattice(std::string const& Description, UnitCell const& uc);

      FiniteLattice(FiniteLattice const& Other);
      FiniteLattice(FiniteLattice&& Other);

      FiniteLattice& operator=(FiniteLattice const& Other);
      FiniteLattice& operator=(FiniteLattice&& Other);

      ~FiniteLattice();

      // returns the number of unit cells in the lattice
      int size() const;

      // returns the total number of sites in the lattice
      int sites() const;

      // returns the n'th unit cell (zero-based)
      UnitCell const& operator[](int i) const;

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

      BraidGroup GetBraidGroup() const;

      // operators

      const_operator_iterator begin_operator() const
      { return Operators_.begin(); }
      const_operator_iterator end_operator() const
      { return Operators_.end(); }

      // returns true if the given operator exits
      bool operator_exists(std::string const& s) const;

      // Sets the description of the operators according to Desc.  Prints a warning to cerr
      // if there are extra operators named in Desc that are not part of the lattice,
      // or if any operators don't have a name.
      void set_operator_descriptions(OperatorDescriptions const& Desc);

      // Lookup the named operator
      FiniteLatticeMPO& operator[](std::string const& Op);
      FiniteLatticeMPO const& operator[](std::string const& Op) const;

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
      FiniteLatticeMPO
      eval_function(Function::OperatorFunction const& Func,
                    Function::ParameterList const& Params) const;

      FiniteLatticeMPO
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

   friend PStream::opstream& operator<<(PStream::opstream& out, FiniteLattice const& L);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, FiniteLattice& L);
};

#endif
