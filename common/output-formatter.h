// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// common/output-formatter.h
//
// Copyright (C) 2016 Ian McCulloch <ian@qusim.net>
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

//
// Basic usage is printf-like,
//
// % [ flag ] [ width ] [ .precision ] [ conversion ]
//
// Allowed conversion characters are d,e,f
// Terminal colours are also allowed, in the form of an escape code
// \[spec]{text}
// where spec is a comma-separated list (may be empty) of either
// integer codes or string descriptive names (not case sensitive)
// as listed in terminal::color.
// A color spec can start with either '\[' (escape code), or with '\\' '[' sequence

#if !defined(MPTOOLKIT_COMMON_OUTPUT_FORMATTER_H)
#define MPTOOLKIT_COMMON_OUTPUT_FORMATTER_H

#include <stdlib.h>
#include <cstring>
#include <string>
#include <algorithm>
#include <cstdio>
#include "common/terminal.h"

enum class ShowColor { Never, Auto, Always };

class OutputFormatter
{
   public:
      OutputFormatter() : DefaultWidth(18), DefaultPrecision(14), Color(false) { }

      explicit OutputFormatter(int DefaultWidth_ = 16,
                               int DefaultPrecision_ = 14,
                               bool Color_ = false) : Color(Color_) { }

      void ShowColors(bool Color_) { Color = Color_; }

      std::string printf(char const* FormatSpec) const;

      template <typename... Args>
      std::string printf(char const* FormatSpec, Args... a) const;

      std::string printf(std::string FormatSpec) const;

      template <typename... Args>
      std::string printf(std::string FormatSpec, Args... a) const;

   private:
      template <typename T, typename... Args>
      void FormatImpl(std::string& Result, char const* p, char const* end,
                      bool InColor, T x, Args... Remain) const;

      template <typename T>
      void FormatImpl(std::string& Result, char const* p, char const* end,
                      bool InColor, T x) const;

      void FormatImpl(std::string& Result, char const* p, char const* end,
                      bool InColor) const;

      int DefaultWidth;
      int DefaultPrecision;
      bool Color;
};

struct FormatSpec
{
   char Flag;  // or \0 if no flag
   int Width;  // 0 if no width specified
   int Precision; // -1 if no precision specified
   char Conversion;  // 'd', 'e', 'f', 'g', 's' for standard printf

   FormatSpec() : Flag('\0'), Width(0), Precision(-1), Conversion('\0') {}

   // precondition: *beg == '%' && *(beg+1) != '%'
   // Parses a FormatSpec at the current location in the string, and
   // returns the FormatSpec object, advances the beg iterator to one-past-the-end
   // of the FormatSpec
   static FormatSpec parse(char const*& beg, char const* end);

   std::string to_string() const;

   template <typename T>
   std::string printf(T x) const;

   static bool is_flag(char c);

};

template <typename T>
std::string
FormatSpec::printf(T x) const
{
   int const bufsz = 1024;
   char buf[bufsz];
   switch (Conversion)
   {
   case 'e' :
   case 'f' :
   case 'g' :
      std::snprintf(buf, bufsz, this->to_string().c_str(), double(x));
      break;
   case 'd' :
      std::snprintf(buf, bufsz, this->to_string().c_str(), int(x));
      break;
   case 's' :
      snprintf(buf, bufsz, this->to_string().c_str(), boost::lexical_cast<std::string>(x).c_str());
      break;
   }
   return std::string(buf);
}

bool
FormatSpec::is_flag(char c)
{
   return c == '0' || c == '+' || c == '#' || c == ' ' || c == '-';
}

FormatSpec
FormatSpec::parse(char const*& beg, char const* end)
{
   FormatSpec Result;
   CHECK_EQUAL(*beg, '%');
   ++beg;
   CHECK(beg < end);
   if (is_flag(*beg))
   {
      Result.Flag = *beg++;
   }
   CHECK(beg < end);
   if (*beg >= '0' && *beg <= '9')
   {
      char* EndPtr;
      Result.Width = std::strtol(beg, &EndPtr, 10);
      beg = EndPtr;
   }
   CHECK(beg < end);
   if (*beg == '.')
   {
      ++beg;
      char* EndPtr;
      Result.Precision = std::strtol(beg, &EndPtr, 10);
      beg = EndPtr;
   }
   CHECK(beg < end);
   Result.Conversion = *beg++;

   TRACE(Result.to_string());

   return Result;
}

std::string
FormatSpec::to_string() const
{
   std::string Out = "%";
   if (Flag != '\0')
      Out += Flag;
   if (Width != 0)
      Out += std::to_string(Width);
   if (Precision != -1)
      Out += '.' + std::to_string(Precision);
   Out += Conversion;
   return Out;
}

template <typename... Args>
inline
std::string
OutputFormatter::printf(std::string Format, Args... x) const
{
   std::string Result;
   FormatImpl(Result, Format.c_str(),
              Format.c_str()+Format.size(), false,
              x...);
   return Result;
}

inline
std::string
OutputFormatter::printf(std::string Format)  const
{
   std::string Result;
   this->FormatImpl(Result, Format.c_str(),
                    Format.c_str()+Format.size(), false);
   return Result;
}

template <typename... Args>
inline
std::string
OutputFormatter::printf(char const* Format, Args... x) const
{
   std::string Result;
   this->FormatImpl(Result, Format,
                    Format + std::strlen(Format), false,
                    x...);
   return Result;
}

inline
std::string
OutputFormatter::printf(char const* Format) const
{
   std::string Result;
   FormatImpl(Result, Format,
              Format+std::strlen(Format), false);
   return Result;
}

template <typename T, typename... Args>
void
OutputFormatter::FormatImpl(std::string& Result,
                            char const* p, char const* end,
                            bool InColor,
                            T x, Args... Remain) const
{
   // text of a \[color]{text} region
   while (*p != '\0')
   {
      if (p[0] == '\[' || (p[0] == '\\' && p[1] == '['))
      {
         // a color spec
         // skip over the initial tokens
         if (p[0] == '\[')
            ++p;
         else
            p += 2;
         char const* e = std::find(p, end, ']');
         std::string ColorSpec = std::string(p,e);
         p = end;
         ++p;
         if (*p != '{')
         {
            abort();
         }
         ++p;
         InColor = true;
         if (Color)
            Result += terminal::parse_color_codes(ColorSpec);
      }
      else if (InColor && p[0] == '\\' && p[1] == '}')
      {
         p += 2;
         Result +=  '}';
      }
      else if (InColor && p[0] == '}')
      {
         ++p;
         if (Color)
            Result += terminal::color_code(terminal::color::Reset);
      }
      else if (p[0] != '%')
      {
         Result += *p++;
      }
      else // p[0] == '%'
      {
         FormatSpec f = FormatSpec::parse(p, end);
         Result += f.printf(x);
         FormatImpl(Result, p, end, InColor, Remain...);
         return;
      }
   }
   // too many parameters supplied
   abort();
}

template <typename T>
void
OutputFormatter::FormatImpl(std::string& Result,
                            char const* p, char const* end,
                            bool InColor,
                            T x) const
{
   // text of a \[color]{text} region
   while (*p != '\0')
   {
      if (p[0] == '\[' || (p[0] == '\\' && p[1] == '['))
      {
         // a color spec
         // skip over the initial tokens
         if (p[0] == '\[')
            ++p;
         else
            p += 2;
         char const* e = std::find(p, end, ']');
         std::string ColorSpec = std::string(p,e);
         p = e;
         ++p;
         if (*p != '{')
         {
            abort();
         }
         ++p;
         InColor = true;
         if (Color)
            Result += terminal::parse_color_codes(ColorSpec);
      }
      else if (InColor && p[0] == '\\' && p[1] == '}')
      {
         p += 2;
         Result +=  '}';
      }
      else if (InColor && p[0] == '}')
      {
         ++p;
         if (Color)
            Result += terminal::color_code(terminal::color::Reset);
      }
      else if (p[0] != '%')
      {
         Result += *p++;
      }
      else // p[0] == '%'
      {
         FormatSpec f = FormatSpec::parse(p, end);
         Result += f.printf(x);
         FormatImpl(Result, p, end, InColor);
         return;
      }
   }
   // too many parameters supplied
   abort();
}

void
OutputFormatter::FormatImpl(std::string& Result,
                            char const* p, char const* end,
                            bool InColor) const
{
   // text of a \[color]{text} region
   while (*p != '\0')
   {
      if (p[0] == '\[' || (p[0] == '\\' && p[1] == '['))
      {
         // a color spec
         // skip over the initial tokens
         if (p[0] == '\[')
            ++p;
         else
            p += 2;
         char const* e = std::find(p, end, ']');
         std::string ColorSpec = std::string(p,e);
         p = end;
         ++p;
         if (*p != '{')
         {
            abort();
         }
         ++p;
         InColor = true;
         if (Color)
            Result += terminal::parse_color_codes(ColorSpec);
      }
      else if (InColor && p[0] == '\\' && p[1] == '}')
      {
         p += 2;
         Result +=  '}';
      }
      else if (InColor && p[0] == '}')
      {
         ++p;
         if (Color)
            Result += terminal::color_code(terminal::color::Reset);
      }
      else if (p[0] != '%')
      {
         Result += *p++;
      }
      else // p[0] == '%'
      {
         abort();
         // too many format specifiers
      }
   }
}

#endif
