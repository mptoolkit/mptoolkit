/* -*- C++ -*- $Id$

  Created 2004-09-05 Ian McCulloch

  originally a 'site block' in DMRG terminology, now a LatticeSite is
  a collection of SiteOperator's defined on some fixed Hilbert space.

  This also has a description, which is just used for informational
  purposes, eg "Fermion site", "Boson site", etc
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

      LatticeSite() : pImpl(new ImplType()) {}

      explicit LatticeSite(std::string const& Description) 
	 : pImpl(new ImplType(Description)) {}

      // precondition: !empty()
      SymmetryList GetSymmetryList() const;

      // precondition: !empty()
      basis1_type const& Basis1() const;
      basis2_type const& Basis2() const;

      bool empty() const { return pImpl->Data.empty(); }

      const_iterator begin() const { return pImpl->Data.begin(); }
      const_iterator end() const { return pImpl->Data.end(); }

      SiteOperator& operator[](std::string const& s) { return pImpl.mutate()->Data[s]; }
      SiteOperator const& operator[](std::string const& s) const;

      const_iterator find(std::string const& s) const { return pImpl->Data.find(s); }

      bool operator_exists(std::string const& s) const { return pImpl->Data.find(s) != pImpl->Data.end(); }

      void CoerceSymmetryList(QuantumNumbers::SymmetryList const& sl);

      std::string const& Description() const { return pImpl->Description; }
      void SetDescription(std::string const& s) { pImpl.mutate()->Description = s; }

      struct ImplType
      {
         std::string Description;
         DataType Data;

	 ImplType() {}
	 ImplType(std::string const& Desc_) : Description(Desc_) {}

         friend PStream::opstream& operator<<(PStream::opstream& out, ImplType const& Impl);
         friend PStream::ipstream& operator>>(PStream::ipstream& in, ImplType& Impl);
      };

   private:
      typedef pvalue_ptr<ImplType> ptr_type;
      ptr_type pImpl;

   friend PStream::opstream& operator<<(PStream::opstream& out, LatticeSite const& B);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, LatticeSite& B);
};

inline
LatticeSite CoerceSL(SymmetryList const& sl, LatticeSite const& s)
{
   LatticeSite r(s);
   r.CoerceSymmetryList(sl);
   return r;
}

// This is used by UnitCell and UnitCellMPO
typedef std::vector<LatticeSite> SiteListType;
typedef pvalue_ptr<SiteListType> SiteListPtrType;
#endif
