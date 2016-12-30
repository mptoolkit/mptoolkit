// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/convertstring.cc
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

#if defined(CONVERTSTRING_NOERRORCHECK)
#include "convertstring_fast.cc"
#else

// conversions to builtins

//
//   The reason why we have two functions here, convert_positive_string and
// convert_negative_string, is to support 2's compliment representations where
// the magnitude of the largest negative number is larger than the largest positive number.
// In that case, attempting a conversion as if it were a positive number would
// produce an overflow.
//   An alternative would be to convert to the corresponding unsigned type, and fix the sign later.
// But then detecting overflow is messier, it needs some traits to get the corresponding unsigned type,
// and the standard might not guarantee that it works.
//
// Two separate functions for +ve and -ve numbers should work properly in all cases.

template <typename T, typename FwdIter>
T
convert_positive_string(FwdIter start, FwdIter end)
{
   if (start == end) throw invalid_string_conversion();
   T n = 0;
   while (start != end)
   {
      typename std::iterator_traits<FwdIter>::value_type c = *start;
      if (c < '0' || c > '9') throw string_conversion_invalid_char();
      T newN = n * 10 + (c - '0');
      if (newN < n) throw string_conversion_overflow();
      n = newN;
      ++start;
   }
   return n;
}

template <typename T, typename FwdIter>
T
convert_negative_string(FwdIter start, FwdIter end)
{
   if (start == end) throw invalid_string_conversion();
   T n = 0;
   while (start != end)
   {
      typename std::iterator_traits<FwdIter>::value_type c = *start;
      if (c < '0' || c > '9') throw string_conversion_invalid_char();
      T newN = n * 10 - (c - '0');
      if (newN > n) throw string_conversion_overflow();
      n = newN;
      ++start;
   }
   return n;
}

// helper function, the default case is is_signed == false
template <typename T, typename FwdIter, bool is_signed = std::numeric_limits<T>::is_signed>
struct convert_string_helper
{
   inline static
   T apply(FwdIter start, FwdIter end)
   {
      return convert_positive_string<T>(start, end);
   }
};

template <typename T, typename FwdIter>
struct convert_string_helper<T, FwdIter, true>
{
   static inline
   T apply(FwdIter start, FwdIter end)
   {
      if (start == end) throw invalid_string_conversion();
      if (*start == '-')
      {
         return convert_negative_string<T>(++start, end);
      }
      // else

      return convert_positive_string<T>(start, end);
   }
};

// convert_string specializations for the builtin types

#define IMPLEMENT_CONVERT_STRING_BUILTIN(type)                          \
template <typename FwdIter>                                             \
struct convert_string_partial<type, FwdIter>                            \
{                                                                       \
   inline static                                                        \
   type                                                                 \
   apply(FwdIter start, FwdIter end)                                    \
   {                                                                    \
     return convert_string_helper<type, FwdIter>::apply(start, end);    \
   }                                                                    \
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

#endif
