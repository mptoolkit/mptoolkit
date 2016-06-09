// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/conflist.h
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

#if !defined(CONFLIST_H_FDG437634YTYJRGE734675TY7REGJ3TY73ER)
#define CONFLIST_H_FDG437634YTYJRGE734675TY7REGJ3TY73ER

#include <map>
#include <string>

#include <sstream>

class ConfList
{
   public:
      // constructs a ConfList.  Any strings of the form ${ENVSTRING} are assumed to be
      // environment strings and are expanded.  If there is no such string, it is left as is.
      // if file does not exist, or cannot be opened for reading, then a std::runtime_error exception
      // is thrown.
      ConfList(std::string const& file);

      template <class T> 
      T Get(const std::string& s, T const& Default) const;

      // overloaded special cases.
      // the std::string version would work with the default, but the special case is much simpler.
      // these should be specializations rather than overloads, but cxx doesn't like the specs here
      // template <> std::string GetConfString<std::string>
      std::string Get(const std::string& s, std::string const& Default) const;

      // This is so it works with string literals as defaults, without specifying std::string
      // as a template parameter.  This is an overload, not a specialization.
      std::string Get(const std::string& s, char const* Default) const;

      // bool is a special case
      // template <> bool GetConfString<bool>
      bool Get(const std::string& s, bool Default) const;

      unsigned long GetBytes(std::string const& s, unsigned long Default) const;

   private:
      std::map<std::string, std::string> Data;
};

template <class T>
T ConfList::Get(const std::string& s, const T& Default) const
{
   std::map<std::string, std::string>::const_iterator I = Data.find(s);
   if (I == Data.end()) return Default;

   std::string Value = I->second;

   std::istringstream InStream(I->second);
   T Val(Default);
 
   InStream >> Val;

   return InStream ? Val : Default;
}

#endif
