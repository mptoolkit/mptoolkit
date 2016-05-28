// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/convertstring.h
//
// Copyright (C) 2000-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

/*
  convertstring.h

  Created 2000-09-01 Ian McCulloch

  Defines the convert_string template function that converts a iterator sequence (of char)
  to some type.  This is supposed to be _fast_, ie there is no default version that
  uses an istringstream or anything like that.  We define conversions for basic
  types here.  Defined types that want a convert_string function should include this header
  and define a specialization of convert_string_partial<T, FwdIter> in their own header.

  the string passed to convert_string must be exactly the right length, with no extraneous whitespace.

  2001-03-29: Its not fast enough.  changed to only do overflow checking if NDEBUG is defined.
  2003-10-20: Reinserted overflow checks, unless CONVERTSTRING_NOERRORCHECK is defined.
              Should work correctly for converting the most negative possible number of signed types.
	      Added overloads for long long.
*/

#if !defined(CONVERTSTRING_H_FHJK348U834UFHJIU45Y789YRUIY34EWI)
#define CONVERTSTRING_H_FHJK348U834UFHJIU45Y789YRUIY34EWI

#include <limits>
#include <stdexcept>
#include <string>
#include <cstring>

//
// partial specializations of convert_string are handled by a static class function.
// By default, the convert_string<T, FwdIter>() function calls convert_string_partial<T, FwdIter>::apply().
// Conversions that encounter invalid chars, too many chars, not enough chars, overflow etc 
// should throw something derived from invalid_string_conversion().
//

class invalid_string_conversion : public std::exception
{
   public:
      invalid_string_conversion()  {}
      char const* what() const throw () { return "Invalid string conversion."; }
};

// exception class for numeric overflow while attempting to convert a string
class string_conversion_overflow : public invalid_string_conversion
{
};

// exception class for too many characters encountered
class string_conversion_too_many_chars : public invalid_string_conversion
{
};

// exception class for not enough characters encountered
class string_conversion_not_enough_chars : public invalid_string_conversion
{
};

// exception class for invalid character encountered
class string_conversion_invalid_char : public invalid_string_conversion
{
};

template <class T, class FwdIter>
struct convert_string_partial
{
   static T apply(FwdIter start, FwdIter end);
};

template <class T, class FwdIter>
inline
T
convert_string(FwdIter start, FwdIter end)
{
   return convert_string_partial<T, FwdIter>::apply(start, end);
}

// a version of convert_string that applies to a std::string

template <class T>
inline
T
convert_string(std::string const& s)
{
   return convert_string<T>(s.begin(), s.end());
}

template <class T>
inline
T
convert_string(char const* s)
{
   return convert_string<T>(s, s+std::strlen(s));
}
#include "convertstring.cc"

#endif
