// -*- C++ -*- $Id$

#include "latticesite.h"

void
CoerceSymmetryList(LatticeSite& site, QuantumNumbers::SymmetryList const& sl)
{
   site.CoerceSymmetryList(sl);
}

void
LatticeSite::CoerceSymmetryList(QuantumNumbers::SymmetryList const& sl)
{
   using ::CoerceSymmetryList;
   ptr_type::lock_type Lock(pImpl.lock());
   for (iterator I = Lock->Data.begin(); I != Lock->Data.end(); ++I)
   {
      CoerceSymmetryList(I->second, sl);
   }
}

SiteOperator const&
LatticeSite::operator[](std::string const& s) const
{ 
   DataType::const_iterator I = pImpl->Data.find(s); 
   CHECK(I != pImpl->Data.end()) << "The site does not contain any operator named " << s;
   return I->second; 
}

SymmetryList
LatticeSite::GetSymmetryList() const
{
   CHECK(!pImpl->Data.empty());
   // if we're debugging, verify that the symmetry list is the same for all operators
#if !defined(NDEBUG)
   SymmetryList ThisSymmetryList = pImpl->Data.begin()->second.GetSymmetryList();
   for (const_iterator I = pImpl->Data.begin(); I != pImpl->Data.end(); ++I)
   {
      CHECK_EQUAL(I->second.GetSymmetryList(), ThisSymmetryList)(I->first);
   }
#endif
   return pImpl->Data.begin()->second.GetSymmetryList();
}

LatticeSite::basis1_type const&
LatticeSite::Basis1() const
{
   CHECK(!pImpl->Data.empty());
   // if we're debugging, verify that the basis is the same for all operators
#if !defined(NDEBUG)
   basis1_type ThisBasis1 = pImpl->Data.begin()->second.Basis1();
   for (const_iterator I = pImpl->Data.begin(); I != pImpl->Data.end(); ++I)
   {
      CHECK_EQUAL(I->second.Basis1(), ThisBasis1);
   }
#endif
   return pImpl->Data.begin()->second.Basis1();
}

LatticeSite::basis2_type const&
LatticeSite::Basis2() const
{
   CHECK(!pImpl->Data.empty());
   // if we're debugging, verify that the basis is the same for all operators
#if !defined(NDEBUG)
   basis2_type ThisBasis2 = pImpl->Data.begin()->second.Basis2();
   for (const_iterator I = pImpl->Data.begin(); I != pImpl->Data.end(); ++I)
   {
      CHECK_EQUAL(I->second.Basis2(), ThisBasis2);
   }
#endif
   return pImpl->Data.begin()->second.Basis2();
}

PStream::opstream& operator<<(PStream::opstream& out, LatticeSite const& B)
{
   return out << B.pImpl;
}

PStream::ipstream& operator>>(PStream::ipstream& in, LatticeSite& B)
{
   return in >> B.pImpl;
}


PStream::opstream& operator<<(PStream::opstream& out, LatticeSite::ImplType const& Impl)
{
   out << Impl.Description;
   out << Impl.Data;
   return out;
}

PStream::ipstream& operator>>(PStream::ipstream& in, LatticeSite::ImplType& Impl)
{
   in >> Impl.Description;
   in >> Impl.Data;
   return in;
}
