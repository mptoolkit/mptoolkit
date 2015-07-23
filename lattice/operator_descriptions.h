// -*- C++ -*- $Id$
//
// Class to hold operator descriptions.
//
// To use:
// Similar to the program_options pattern,
// OperatorDescrptions OpDescriptions;
// OpDescriptions.add_operators()("First operator", "First operator description")("second operator", "description) ... ;

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

   private:
      data_type Descriptions;
      std::map<std::string, int> Index;

};

inline
std::ostream& operator<<(std::ostream& out, OperatorDescriptions const& d)
{
   for (OperatorDescriptions::const_iterator I = d.begin(); I != d.end(); ++I)
   {
      out << std::setw(10) << std::left << I->first << " - " << I->second << '\n';
   }
   return out;
}

#endif
