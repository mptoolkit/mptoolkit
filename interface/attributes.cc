// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// interface/attributes.cc
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
// Attribute
//

template <typename T>
Attribute::Attribute(T const& x)
   : value_(boost::lexical_cast<std::string>(x))
{
}

template <typename T>
Attribute& Attribute::operator=(T const& x)
{
   value_ = boost::lexical_cast<std::string>(x);
   return *this;
}

template <typename T>
Attribute& Attribute::operator+=(T const& x)
{
   value_ += boost::lexical_cast<std::string>(x);
   return *this;
}

#if 0
template <typename T>
Attribute::operator T() const
{
   return boost::lexical_cast<T>(value_);
}
#endif

template <typename T>
T Attribute::as() const
{
   return boost::lexical_cast<T>(value_);
}

template <typename T>
T Attribute::get_or_default(T const& Def) const
{
   if (value_.empty())
      return Def;
   return boost::lexical_cast<T>(value_);
}
