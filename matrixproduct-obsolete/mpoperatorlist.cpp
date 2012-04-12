// -*- C++ -*- $Id$

#include "mpoperatorlist.h"

OperatorList::OperatorList(Lattice const& L_) : L(L_)
{
#if 0
   // this is now redundant, we construct site operators on the fly
   for (int i = 0; i < L.size(); ++i)
   {
      std::string Coords = '(' + L.Coords(i) + ')';
      for (SiteBlock::const_iterator I = L.BlockAtSite(i).begin();
	   I != L.BlockAtSite(i).end(); ++I)
      {
	 if (I->first == "I") continue;  // don't bother adding the identity operator!
	 Data[I->first + Coords] = CreateMPOperator(L, I->first, i);
      }
   }
#endif
   // now we add (just one!) identity operator for the whole lattice
   Data["I"] = CreateMPOperator(L, "I", 0);
}

OperatorList::OperatorType& 
OperatorList::operator[](std::string const& s) 
{ 
   if (Data.find(s) == Data.end())
      Data[s] = this->DoConstruct(s);
   return Data[s]; 
}

OperatorList::OperatorType
OperatorList::DoConstruct(std::string const& Operator) const
{
   std::string::const_iterator I = std::find(Operator.begin(), Operator.end(), '(');
   if (I == Operator.end())
      return OperatorType();
   std::string OpName(Operator.begin(), I);
   ++I;

   std::string::const_iterator J = Operator.end(); --J;
   if (*J != ')')
      return OperatorType();
   std::string sSite(I,J);

   int nSite = this->GetLattice().site_at_coordinate(sSite);
   if (nSite < 0)
   {
      return OperatorType();
   }

   Lattice::const_iterator K = this->GetLattice().begin();
   std::advance(K, nSite-1);
   if (K->find(OpName) == K->end())
      return OperatorType();

   return CreateMPOperator(this->GetLattice(), OpName, nSite-1);
}

OperatorList::OperatorType
OperatorList::operator[](std::string const& Operator) const
{
   const_iterator dI = Data.find(Operator);
   if (dI != Data.end())
      return dI->second;

   std::string::const_iterator I = std::find(Operator.begin(), Operator.end(), '(');
   if (I == Operator.end())
   {
      PANIC("Operator not found")(Operator);
   }
   std::string OpName(Operator.begin(), I);
   ++I;

   std::string::const_iterator J = Operator.end(); --J;
   if (*J != ')')
   {
      PANIC("Operator not found")(Operator);
   }
   std::string sSite(I,J);

   int nSite = this->GetLattice().site_at_coordinate(sSite);
   if (nSite < 0)
   {
      PANIC("Operator not found")(Operator)(sSite)(OpName);
   }

   return CreateMPOperator(this->GetLattice(), OpName, nSite-1);
}

bool OperatorList::HasOperator(std::string const& Operator) const
{
   if (Data.find(Operator) != Data.end())
      return true;

   std::string::const_iterator I = find(Operator.begin(), Operator.end(), '(');
   if (I == Operator.end())
      return false;
   std::string OpName(Operator.begin(), I);
   ++I;

   std::string::const_iterator J = Operator.end(); --J;
   if (*J != ')')
      return false;
   std::string sSite(I,J);

   int nSite = this->GetLattice().site_at_coordinate(sSite);
   if (nSite < 0)
      return false;

   Lattice::const_iterator K = this->GetLattice().begin();
   std::advance(K, nSite-1);
   return (K->find(OpName) != K->end());
}

PStream::opstream& operator<<(PStream::opstream& out, OperatorList const& L)
{
   return out << L.Data << L.L;
}

PStream::ipstream& operator>>(PStream::ipstream& in, OperatorList& L)
{
   return in >> L.Data >> L.L;
}

