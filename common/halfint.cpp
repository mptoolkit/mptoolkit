// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/halfint.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "halfint.h"
#include <iomanip>
#include <sstream>

std::ostream& operator<<(std::ostream& out, const half_int& H)
{
   out << (H.twice() / 2.0);
   return out;
}

std::istream& operator>>(std::istream& in, half_int& H)
{
#if 0
   // this is a (possibly) faster (except for the string allocation), but less robust version
   in >> skipws;
   char c;
   std::string s;
   in >> c;
   s = c;
   TRACE(c);
   while ((in >> c) && ((c >= '0' && c <= '9') || c == '.' || c == '/'))
   {
      s += c;
   }
   in.putback(c);
   TRACE(s);
   H = convert_string<half_int>(s.begin(), s.end());
   return in;
#else
   double d;
   in >> d;
   H = half_int(d);
   return in;
#endif
}

void half_int::throw_cannot_convert()
{
   throw std::runtime_error("half_int: cannot convert odd half_int to integer!");
}

std::string to_string_fraction(half_int h)
{
   std::ostringstream s;
   if (is_integral(h))
   {
      s << to_int_assert(h);
   }
   else
   {
      s << h.twice() << "/2";
   }
   return s.str();
}
