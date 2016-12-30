// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/environment.h
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

#if !defined(ENVIRONMENT_H_KHCJKSHDUY478Y5478YOPEW)
#define ENVIRONMENT_H_KHCJKSHDUY478Y5478YOPEW

#include <string>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <stdlib.h> // for getenv

inline
std::string quote_shell(std::string s)
{
   if (!boost::all(s, !boost::is_any_of(" \t()[]*\\\"")))
   {
      boost::replace_all(s, std::string("\\"), std::string("\\\\"));
      boost::replace_all(s, std::string("\""), std::string("\\\""));
      return '"' + s + '"';
   }
   // else
   return s;
}

inline
std::string cmdline(int argc, char** argv)
{
   std::string result;
   if (argc > 0) result = quote_shell(argv[0]);
   for (int i = 1; i < argc; ++i)
      result = result + ' ' + quote_shell(argv[i]);
   return result;
}

// expands components of the form ${XXX} as the corresponding environment string
std::string ExpandEnvironment(std::string const& s);

// getenv_or_default
// Get an environment string, or if the environment string is not defined, return the specified Default value
template <typename T>
T getenv_or_default(std::string const& str, T const& Default);

template <typename T>
T getenv_or_default(std::string const& str, T const& Default)
{
   try
   {
      char const* Str = getenv(str.c_str());
      if (Str != NULL)
         return boost::lexical_cast<T>(Str);
   }
   catch (boost::bad_lexical_cast&)
   {
   }
   catch (...)
   {
      throw;
   }

   return Default;
}

#endif
