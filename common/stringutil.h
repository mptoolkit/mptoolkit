// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/stringutil.h
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

/*
  stringutil.h

  Various utilities for manipulating strings

  First assembled from various fragments 2000-07-22 Ian McCulloch
*/

#if !defined(STRINGUTIL_H_SFDHJ38UFJ38UCNUIFJIC)
#define STRINGUTIL_H_SFDHJ38UFJ38UCNUIFJIC

#include <stdexcept>
#include <string>
#include <ctype.h>
#include <sstream>
#include <algorithm>
#include "common/trace.h"

//
// RemoveWhiteSpace removes white space at the front and the back of the string, according to isspace()
//
inline void RemoveWhiteSpace(std::string& s)
{
   while (s.length() > 1 && isspace(s[0])) s.erase(0, 1);
   while (s.length() > 1 && isspace(s[s.length()-1])) s.erase(s.length()-1, 1);
}

//
// RemoveAllWhiteSpace removes all white space characters, even internal white space.
//

inline
void RemoveAllWhiteSpace(std::string& s)
{
   size_t i = 0;
   while (i < s.size())
   {
      if (isspace(s[i])) s.erase(i, 1);
      else ++i;
   }
}

//
// ConvertString is a template function that converts a string to an arbitary type via a istringstream.
//

template <class T>
T
ConvertString(const char* str)
{
   CHECK(str != NULL);
   std::istringstream Stream(str);
   Stream.exceptions(std::ios_base::badbit|std::ios_base::failbit);
   T Temp;
   try
   {
      Stream >> Temp;
   }
   catch (std::ios_base::failure& f)
   {
      throw std::runtime_error(std::string(f.what()) + " converting " + str + " to " + typeid(T).name());
   }
   return Temp;
}

template <class T>
T
ConvertString(std::string const& str)
{
   std::istringstream Stream(str);
   Stream.exceptions(std::ios_base::badbit|std::ios_base::failbit);
   T Temp;
   try
   {
      Stream >> Temp;
   }
   catch (std::ios_base::failure& f)
   {
      throw std::runtime_error(std::string(f.what()) + " converting " + str + " to " + typeid(T).name());
   }
   return Temp;
}

//
// ConvertToString is the inverse of ConvertString - conversion of some type to a string via a ostringstream
//

template <class T>
std::string
ConvertToString(T const& x)
{
   std::ostringstream Stream;
   Stream << x;
   return Stream.str();
}

// splits a character string into substrings separated by Delim.
template <class VType, class OutIter2>
void
Split(std::string const& q, VType Delim, OutIter2 OutIter)
{
   // [q_beg, q_end) is a subsequence terminated by Delim
   std::string::const_iterator q_beg = q.begin();
   std::string::const_iterator q_end = std::find(q_beg, q.end(), Delim);

   // loop until our subsequence is terminated by q.end()
   while (q_end != q.end())
   {
      *OutIter = std::string(q_beg, q_end);
      ++OutIter;

      q_beg = q_end;
      ++q_beg;          // skip over Delim

      q_end = std::find(q_beg, q.end(), Delim);
   }

   // output the final subsequence
   if (q_beg != q_end)
   {
      *OutIter = std::string(q_beg, q_end);
      ++OutIter;
   }
}


#endif
