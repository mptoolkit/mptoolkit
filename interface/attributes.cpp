// -*- C++ -*- $Id$

#include "attributes.h"

//
// Attribute
//

Attribute::Attribute()
   : value_()
{
}

std::ostream& operator<<(std::ostream& out, Attribute const& a)
{
   return out << a.value_;
}

PStream::opstream& operator<<(PStream::opstream& out, Attribute const& a)
{
   //   using PStream::operator<<;
   return out << a.value_;
}

PStream::ipstream& operator>>(PStream::ipstream& in, Attribute& a)
{
   //   using PStream::operator>>;
   return in >> a.value_;
}

//
// AttributeList
//

std::ostream& operator<<(std::ostream& out, AttributeList const& l)
{
   for (AttributeList::const_iterator I = l.begin(); I != l.end(); ++I)
      out << I->first << '=' << I->second << '\n';
   return out;
}

PStream::opstream& operator<<(PStream::opstream& out, AttributeList const& l)
{
   return out << l.data_;
}

PStream::ipstream& operator>>(PStream::ipstream& in, AttributeList& l)
{
   return in >> l.data_;
}

Attribute 
AttributeList::operator[](std::string const& s) const
{
   const_iterator I = data_.find(s);
   if (I != data_.end()) return I->second;
   // else
   return Attribute();
}

