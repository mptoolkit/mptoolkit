// -*- C++ -*- $Id$

#if !defined(MPOPERATORLIST_H_JHDSCIUREWHFUHFYT87YF87YO)
#define MPOPERATORLIST_H_JHDSCIUREWHFUHFYT87YF87YO

#include "matrixproduct/mpoperator.h"
#include "lattice.h"
#include <boost/lexical_cast.hpp>

class OperatorList
{
   private:
      typedef std::map<std::string, MPOperator> DataType;

   public:
      typedef MPOperator               OperatorType;
      typedef DataType::value_type     value_type;
      typedef DataType::iterator       iterator;
      typedef DataType::const_iterator const_iterator;

      OperatorList() {}

      OperatorList(Lattice const& L);

      int size() const { return L.size(); }

      iterator begin() { return Data.begin(); }
      iterator end() { return Data.end(); }

      const_iterator begin() const { return Data.begin(); }
      const_iterator end() const { return Data.end(); }

      OperatorType& operator[](std::string const& s);
      OperatorType operator[](std::string const& s) const;

      OperatorType Lookup(std::string const& Operator, std::string const& Site) const
	 { return this->operator[](Operator + '(' + Site +')'); }

      bool HasOperator(std::string const& Operator) const;

      bool HasOperator(std::string const& Operator, std::string const& Site) const
	 { return this->HasOperator(Operator + '(' + Site + ')'); }

      Lattice const& GetLattice() const { return L; }

      SymmetryList GetSymmetryList() const { return L.front().GetSymmetryList(); }

   private:
      OperatorType DoConstruct(std::string const& Operator) const;

      DataType Data;
      Lattice L;

   friend PStream::opstream& operator<<(PStream::opstream& out, OperatorList const& L);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, OperatorList& L);
};


#endif
