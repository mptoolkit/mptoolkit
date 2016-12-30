// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// pstream/optional.h
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
/* -*- C++ -*- $Id$

  Seralization of boost::optional
*/

#if !defined(OPTIONAL_H_JSDCHWUI4357784Y7WEHOLWEHO)
#define OPTIONAL_H_JSDCHWUI4357784Y7WEHOLWEHO

#include "pstream.h"
#include <boost/optional.hpp>

namespace PStream
{

template <int format, typename T>
opstreambuf<format>& operator<<(opstreambuf<format>& out, boost::optional<T> const& x)
{
   bool b = x;
   out << b;
   if (b)
      out << (*x);
   return out;
}

template <int format, typename T>
ipstreambuf<format>& operator>>(ipstreambuf<format>& in, boost::optional<T>& x)
{
   bool b; in >> b;
   if (b)
   {
      T n;
      in >> n;
      x = n;
   }
   else
      x = boost::optional<T>();
   return in;
}

} // namespace PStream

#endif
