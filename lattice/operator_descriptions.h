// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// lattice/operator_descriptions.h
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
//
// Class to hold operator descriptions.
//
// To use:
// Similar to the program_options pattern,
// OperatorDescriptions OpDescriptions;
// OpDescriptions.add_operators()("First operator", "First operator description")
//                               ("second operator", "description) ... ;
//
// functions:
// OpDescriptions.func("Name", "Description")(args) = "definition";
//

#if !defined(MPTOOLKIT_LATTICE_OPERATOR_DESCRIPTIONS_H)
#define MPTOOLKIT_LATTICE_OPERATOR_DESCRIPTIONS_H

#include <string>
#include <map>
#include <iostream>
#include <iomanip>

class OperatorDescriptions
{
   public:
      typedef std::pair<std::string, std::string> value_type;
      typedef std::vector<value_type> data_type;

      typedef data_type::const_iterator const_iterator;

      typedef std::map<std::string, std::string> function_list_type;
      typedef function_list_type::const_iterator const_function_iterator;

   private:
      typedef std::map<std::string, int> index_type;

   public:
      typedef std::vector<std::pair<std::string, std::string>> authors_type;

      std::string description() const { return Description; }

      void description(std::string s) { Description = std::move(s); }

      void author(std::string Name, std::string Email) { Authors.push_back({Name, Email}); }

      authors_type const& authors() const { return Authors; }

      unsigned author_size() const { return Authors.size(); }
      std::pair<std::string, std::string> author(int n) const { return Authors[n]; }

      unsigned size() const { return Descriptions.size(); }

      const_iterator begin() const { return Descriptions.begin(); }
      const_iterator end() const { return Descriptions.end(); }

      struct OperatorDescProxy
      {
	 OperatorDescProxy(data_type& Desc_, index_type& Index_)
	    : Descriptions(Desc_), Index(Index_) {}

	 OperatorDescProxy const& operator()(std::string const& Name,
					     std::string const& Desc) const
	 {
	    if (Index.find(Name) == Index.end())
	    {
	       Index[Name] = Descriptions.size();
	       Descriptions.push_back(std::make_pair(Name, Desc));
	    }
	    else
	    {
	       Descriptions[Index[Name]] = std::make_pair(Name, Desc);
	    }
	    return *this;
	 }
	 data_type& Descriptions;
	 index_type& Index;
      };

      OperatorDescProxy add_operators() { return OperatorDescProxy(Descriptions, Index); }

      // unit cell operators

      const_iterator cell_begin() const { return CellOperatorDescriptions.begin(); }
      const_iterator cell_end() const { return CellOperatorDescriptions.end(); }

      OperatorDescProxy add_cell_operators() { return OperatorDescProxy(CellOperatorDescriptions, 
									CellIndex); }

      // Functions

      struct FunctionDescProxy
      {
	 FunctionDescProxy(function_list_type& F_)
	    : F(&F_) {}

	 FunctionDescProxy const& operator()(std::string const& Name,
					     std::string const& Desc) const
	 {
	    (*F)[Name] = Desc;
	    return *this;
	 }

	 function_list_type* F;
      };

      FunctionDescProxy add_functions() { return FunctionDescProxy(Functions); }

      unsigned size_function() const { return Functions.size(); }

      const_function_iterator begin_function() const { return Functions.begin(); }
      const_function_iterator end_function() const { return Functions.end(); }

   private:
      std::string Description;
      std::vector<std::pair<std::string, std::string>> Authors;

      data_type Descriptions;
      std::map<std::string, int> Index;

      data_type CellOperatorDescriptions;
      std::map<std::string, int> CellIndex;

      function_list_type Functions;
};

inline
std::ostream& operator<<(std::ostream& out, OperatorDescriptions const& d)
{
   out << "Description: " << d.description() << '\n';
   for (auto const& x : d.authors())
   {
      out << "Author: " << x.first << " <" << x.second << ">\n";
   }
   out << "Operators:\n";
   if (d.size() == 0)
      out << "(none)\n";
   for (OperatorDescriptions::const_iterator I = d.begin(); I != d.end(); ++I)
   {
      out << std::setw(10) << std::left << I->first << " - " << I->second << '\n';
   }

   out << "\nFunctions:\n";
   if (d.size_function() == 0)
      out << "(none)\n";
   for (OperatorDescriptions::const_function_iterator I = d.begin_function(); 
	I != d.end_function(); ++I)
   {
      out << std::setw(10) << std::left << I->first << " - " << I->second << '\n';
   }

   return out;
}

#endif
