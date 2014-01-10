/* -*- C++ -*- $Id$

  Created 2004-09-05 Ian McCulloch

  originally a 'site block' in DMRG terminology, now a LatticeSite is
  a collection of SiteOperator's defined on some fixed Hilbert space.
*/

#if !defined(LATTICESITE_H_853YF987RHVHYC85HYT87HGGO2)
#define LATTICESITE_H_853YF987RHVHYC85HYT87HGGO2

#include "pstream/pstream.h"
#include "pheap/pvalueptr.h"
#include "quantumnumbers/symmetrylist.h"
#include "siteoperator.h"
#include <map>

using QuantumNumbers::SymmetryList;

class LatticeSite
{
   private:
      typedef std::map<std::string, SiteOperator> DataType;

   public:
      typedef DataType::value_type        value_type;
      typedef DataType::iterator          iterator;
      typedef DataType::const_iterator    const_iterator;
      typedef SiteOperator::basis1_type   basis1_type;
      typedef SiteOperator::basis2_type   basis2_type;

      LatticeSite() : Data(new DataType()) {}

      // precondition: !empty()
      SymmetryList GetSymmetryList() const;

      // precondition: !empty()
      basis1_type const& Basis1() const;
      basis2_type const& Basis2() const;

      bool empty() const { return Data->empty(); }

      const_iterator begin() const { return Data->begin(); }
      const_iterator end() const { return Data->end(); }

      SiteOperator& operator[](std::string const& s) { return (*Data.mutate())[s]; }
      SiteOperator const& operator[](std::string const& s) const;

      const_iterator find(std::string const& s) const { return Data->find(s); }

      bool exists(std::string const& s) const { return Data->find(s) != Data->end(); }

      void CoerceSymmetryList(QuantumNumbers::SymmetryList const& sl);

   private:
      typedef pvalue_ptr<DataType> ptr_type;
      ptr_type Data;

   friend PStream::opstream& operator<<(PStream::opstream& out, LatticeSite const& B);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, LatticeSite& B);
};

#endif
