// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// interface/attributes.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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
