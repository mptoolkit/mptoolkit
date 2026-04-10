// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// common/json_emit.h
//
// Copyright (C) 2026 Ian McCulloch <ian@qusim.net>
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

#if !defined(MPTOOLKIT_COMMON_JSON_EMIT_H)
#define MPTOOLKIT_COMMON_JSON_EMIT_H

#include <ostream>
#include <string>

namespace json_emit
{

inline
void emit_string(std::ostream& Out, std::string const& Str)
{
   Out << '"';
   for (unsigned char ch : Str)
   {
      switch (ch)
      {
      case '\"':
         Out << "\\\"";
         break;
      case '\\':
         Out << "\\\\";
         break;
      case '\b':
         Out << "\\b";
         break;
      case '\f':
         Out << "\\f";
         break;
      case '\n':
         Out << "\\n";
         break;
      case '\r':
         Out << "\\r";
         break;
      case '\t':
         Out << "\\t";
         break;
      default:
         if (ch < 0x20)
         {
            static char const Hex[] = "0123456789abcdef";
            Out << "\\u00" << Hex[(ch >> 4) & 0x0f] << Hex[ch & 0x0f];
         }
         else
         {
            Out << ch;
         }
         break;
      }
   }
   Out << '"';
}

inline
void emit_key(std::ostream& Out, std::string const& Key)
{
   emit_string(Out, Key);
   Out << ':';
}

class object
{
public:
   explicit object(std::ostream& Out_)
      : Out(Out_),
        First(true)
   {
      Out << '{';
   }

   ~object()
   {
      Out << '}';
   }

   void next(std::string const& Key)
   {
      if (!First)
         Out << ',';
      First = false;
      emit_key(Out, Key);
   }

private:
   std::ostream& Out;
   bool First;
};

class array
{
public:
   explicit array(std::ostream& Out_)
      : Out(Out_),
        First(true)
   {
      Out << '[';
   }

   ~array()
   {
      Out << ']';
   }

   void next()
   {
      if (!First)
         Out << ',';
      First = false;
   }

private:
   std::ostream& Out;
   bool First;
};

} // namespace json_emit

#endif
