// -*- C++ -*- $Id$
//
// Class to hold operator descriptions.
//
// To use:
// Similar to the program_options pattern,
// OperatorDescrptions OpDescriptions;
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

      unsigned size() const { return Descriptions.size(); }

      const_iterator begin() const { return Descriptions.begin(); }
      const_iterator end() const { return Descriptions.end(); }

      OperatorDescriptions& operator()(std::string const& Name,
				       std::string const& Desc)
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
      OperatorDescriptions& add_operators() { return *this; }

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
      data_type Descriptions;
      std::map<std::string, int> Index;
      function_list_type Functions;
};

inline
std::ostream& operator<<(std::ostream& out, OperatorDescriptions const& d)
{
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
