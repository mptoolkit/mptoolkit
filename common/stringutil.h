// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// common/stringutil.h
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
#include <cctype>
#include <sstream>
#include <algorithm>
#include <string_view>
#include <vector>
#include "common/trace.h"

//
// RemoveWhiteSpace removes white space at the front and the back of the string, according to isspace()
//
inline void RemoveWhiteSpace(std::string& s)
{
   while (s.length() > 1 && isspace(s[0])) s.erase(0, 1);
   while (s.length() > 1 && isspace(s[s.length()-1])) s.erase(s.length()-1, 1);
}

inline
bool StartsWith(std::string_view s, std::string_view prefix)
{
   return s.starts_with(prefix);
}

inline
bool IEquals(std::string_view x, std::string_view y)
{
   if (x.size() != y.size())
      return false;

   for (std::string::size_type i = 0; i < x.size(); ++i)
   {
      if (std::tolower(static_cast<unsigned char>(x[i])) !=
          std::tolower(static_cast<unsigned char>(y[i])))
         return false;
   }

   return true;
}

inline
void ReplaceAll(std::string& s, std::string_view from, std::string_view to)
{
   if (from.empty())
      return;

   std::string::size_type pos = 0;
   while ((pos = s.find(from.data(), pos, from.size())) != std::string::npos)
   {
      s.replace(pos, from.size(), to.data(), to.size());
      pos += to.size();
   }
}

inline
std::string TrimCopy(std::string_view s)
{
   auto const is_not_space = [](unsigned char c) { return !std::isspace(c); };
   auto first = std::find_if(s.begin(), s.end(), is_not_space);
   auto last = std::find_if(s.rbegin(), s.rend(), is_not_space).base();

   if (first >= last)
      return std::string();

   return std::string(first, last);
}

inline
void Trim(std::string& s)
{
   s = TrimCopy(s);
}

inline
std::string ToLowerCopy(std::string s)
{
   std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
   return s;
}

inline
std::vector<std::string> SplitCompress(std::string_view s, std::string_view separators)
{
   std::vector<std::string> result;
   std::string_view::size_type begin = 0;
   while (begin < s.size())
   {
      begin = s.find_first_not_of(separators, begin);
      if (begin == std::string_view::npos)
         break;

      std::string_view::size_type end = s.find_first_of(separators, begin);
      result.emplace_back(s.substr(begin, end - begin));
      if (end == std::string_view::npos)
         break;
      begin = end + 1;
   }

   return result;
}

template <typename Range>
std::string JoinStrings(Range const& strings, std::string_view separator)
{
   std::ostringstream out;
   for (auto i = strings.begin(); i != strings.end(); ++i)
   {
      if (i != strings.begin())
         out << separator;
      out << *i;
   }
   return out.str();
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

template <>
inline
std::string
ConvertString<std::string>(char const* str)
{
   CHECK(str != NULL);
   return str;
}

template <>
inline
std::string
ConvertString<std::string>(std::string const& str)
{
   return str;
}

//
// ConvertStringStrict is for user-facing or persisted text where malformed
// trailing characters should be rejected.
//

template <class T>
T
ConvertStringStrict(const char* str)
{
   CHECK(str != NULL);
   std::istringstream Stream(str);
   Stream.exceptions(std::ios_base::badbit|std::ios_base::failbit);
   T Temp;
   try
   {
      Stream >> Temp;
      Stream.exceptions(std::ios_base::badbit);
      Stream >> std::ws;
   }
   catch (std::ios_base::failure& f)
   {
      throw std::runtime_error(std::string(f.what()) + " converting " + str + " to " + typeid(T).name());
   }
   if (!Stream.eof())
      throw std::runtime_error("trailing characters converting " + std::string(str) + " to " + typeid(T).name());
   return Temp;
}

template <class T>
T
ConvertStringStrict(std::string const& str)
{
   std::istringstream Stream(str);
   Stream.exceptions(std::ios_base::badbit|std::ios_base::failbit);
   T Temp;
   try
   {
      Stream >> Temp;
      Stream.exceptions(std::ios_base::badbit);
      Stream >> std::ws;
   }
   catch (std::ios_base::failure& f)
   {
      throw std::runtime_error(std::string(f.what()) + " converting " + str + " to " + typeid(T).name());
   }
   if (!Stream.eof())
      throw std::runtime_error("trailing characters converting " + str + " to " + typeid(T).name());
   return Temp;
}

template <>
inline
std::string
ConvertStringStrict<std::string>(char const* str)
{
   CHECK(str != NULL);
   return str;
}

template <>
inline
std::string
ConvertStringStrict<std::string>(std::string const& str)
{
   return str;
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
