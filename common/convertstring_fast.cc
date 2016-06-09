// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/convertstring_fast.cc
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
  convertstring_fast.cc

  fast, no error checking implementation of convertstring.cc, convertstring for builtins.

  This version is depreciated.
*/

template <class T, class FwdIter>
T
convert_positive_string(FwdIter start, FwdIter end)
{
   T n = 0;
   while (start != end)
   {
      n = n * 10 + (*start - '0');
      ++start;
   }
   return n;
}

// helper function, the default case is is_signed == false
template <class T, class FwdIter, bool is_signed = std::numeric_limits<T>::is_signed>
struct convert_string_helper
{
   inline static
   T apply(FwdIter start, FwdIter end)
   {
      return convert_positive_string<T>(start, end);
   }
};

template <class T, class FwdIter>
struct convert_string_helper<T, FwdIter, true>
{
   static inline
   T apply(FwdIter start, FwdIter end)
   {
      if (*start == '-') 
      {
         return -convert_positive_string<T>(++start, end);
      }
      // else

      return convert_positive_string<T>(start, end);
   }
};

// convert_string specializations for the builtin types

#define IMPLEMENT_CONVERT_STRING_BUILTIN(type)				\
template <class FwdIter>						\
struct convert_string_partial<type, FwdIter>				\
{									\
   inline static       							\
   type									\
   apply(FwdIter start, FwdIter end)					\
   {									\
     return convert_string_helper<type, FwdIter>::apply(start, end);	\
   }									\
};

IMPLEMENT_CONVERT_STRING_BUILTIN(char)
IMPLEMENT_CONVERT_STRING_BUILTIN(signed char)
IMPLEMENT_CONVERT_STRING_BUILTIN(unsigned char)
IMPLEMENT_CONVERT_STRING_BUILTIN(short)
IMPLEMENT_CONVERT_STRING_BUILTIN(unsigned short)
IMPLEMENT_CONVERT_STRING_BUILTIN(int)
IMPLEMENT_CONVERT_STRING_BUILTIN(unsigned int)
IMPLEMENT_CONVERT_STRING_BUILTIN(long)
IMPLEMENT_CONVERT_STRING_BUILTIN(unsigned long)
#if defined(USE_LONGLONG)
IMPLEMENT_CONVERT_STRING_BUILTIN(long long)
IMPLEMENT_CONVERT_STRING_BUILTIN(unsigned long long)
#endif
