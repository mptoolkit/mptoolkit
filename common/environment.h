// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// common/environment.h
//
// Copyright (C) 2004-2023 Ian McCulloch <ian@qusim.net>
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

#if !defined(ENVIRONMENT_H_KHCJKSHDUY478Y5478YOPEW)
#define ENVIRONMENT_H_KHCJKSHDUY478Y5478YOPEW

#include "common/stringutil.h"
#include <ctime>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <stdlib.h> // for getenv
#include <sys/stat.h>

inline
std::string quote_shell(std::string s)
{
   if (s.find_first_of(" \t()[]*\\\"") != std::string::npos)
   {
      ReplaceAll(s, std::string("\\"), std::string("\\\\"));
      ReplaceAll(s, std::string("\""), std::string("\\\""));
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

inline
bool env_exists(std::string const& str)
{
   return getenv(str.c_str()) != nullptr;
}

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
         return ConvertStringStrict<T>(Str);
   }
   catch (std::runtime_error const&)
   {
   }
   catch (...)
   {
      throw;
   }

   return Default;
}

inline
char const*
getenv_or_default(std::string const& str, char const* Default)
{
   char const* Str = getenv(str.c_str());
   return Str ? Str : Default;
}

// Format the current local date/time using strftime format specifiers.
inline
std::string format_date_time(std::string const& Format)
{
   std::time_t Time = std::time(nullptr);
   std::tm LocalTime;
#if defined(_WIN32) || defined(_MSC_VER)
   if (localtime_s(&LocalTime, &Time) != 0)
      return std::string();
#else
   if (localtime_r(&Time, &LocalTime) == nullptr)
      return std::string();
#endif

   std::ostringstream Out;
   Out << std::put_time(&LocalTime, Format.c_str());
   return Out ? Out.str() : std::string();
}

inline
bool bench_filename_exists(std::string const& FileName)
{
   struct stat Buffer;
   return stat(FileName.c_str(), &Buffer) == 0;
}

inline
bool parse_bench_filename_width(std::string const& WidthStr, int& Width)
{
   int const MaxWidth = 100;
   Width = 0;
   for (char c : WidthStr)
   {
      if (c < '0' || c > '9')
         return false;
      int const Digit = c - '0';
      // Avoid exceptions and pathological filename padding from bad environment values.
      if (Width > (MaxWidth - Digit) / 10)
         return false;
      Width = Width * 10 + Digit;
   }
   return true;
}

// Find an unused filename of the form BaseName + N + Extension, where N is
// zero-padded to at least Width digits.
inline
std::string find_next_unused_filename(std::string const& BaseName, std::string const& Extension,
                                      int Width = 0)
{
   int Number = 1;
   std::string FileName;
   std::ostringstream Out;
   do
   {
      Out.str("");
      Out.clear();
      Out << std::setw(Width) << std::setfill('0') << Number++;
      FileName = BaseName + Out.str() + Extension;
   } while (bench_filename_exists(FileName));
   return FileName;
}

// Replace filename format specifiers used by MP_BENCHFILE.
inline
std::string replace_filename_format_specifiers(std::string FileName)
{
   struct FormatSpecifier
   {
      char const* Token;
      char const* Format;
   };

   FormatSpecifier const FormatMap[] = {
      {"%D", "%Y-%m-%d-%H:%M:%S"},
      {"%T", "%H:%M:%S"},
      {"%F", "%Y-%m-%d"},
      {"%C", "%Y%m%d_%H%M%S"},
      {"%z", "%z"},
      {"%Z", "%Z"}
   };

   for (auto const& Format : FormatMap)
   {
      std::size_t Pos = 0;
      std::size_t const TokenLength = std::char_traits<char>::length(Format.Token);
      while ((Pos = FileName.find(Format.Token, Pos)) != std::string::npos)
      {
         std::string Replacement = format_date_time(Format.Format);
         FileName.replace(Pos, TokenLength, Replacement);
         Pos += Replacement.length();
      }
   }

   std::size_t Pos = FileName.find('%');
   while (Pos != std::string::npos)
   {
      std::size_t End = Pos + 1;
      while (End < FileName.size() && FileName[End] >= '0' && FileName[End] <= '9')
         ++End;
      if (End < FileName.size() && FileName[End] == 'n')
      {
         std::string WidthStr = FileName.substr(Pos + 1, End - Pos - 1);
         int Width = 0;
         if (!parse_bench_filename_width(WidthStr, Width))
            return std::string();
         FileName = find_next_unused_filename(FileName.substr(0, Pos), FileName.substr(End + 1),
                                              Width);
         break;
      }
      Pos = FileName.find('%', Pos + 1);
   }

   return FileName;
}

// Open an MP_BENCHFILE-style output stream. A leading '+' appends, and a
// leading '++' appends after writing a blank separator line.
inline
std::ofstream open_bench_file(std::string FileName)
{
   if (FileName.empty())
      return std::ofstream();

   bool Append = false;
   bool BlankLine = false;
   if (!FileName.empty() && FileName[0] == '+')
   {
      Append = true;
      FileName.erase(0, 1);
      if (!FileName.empty() && FileName[0] == '+')
      {
         BlankLine = true;
         FileName.erase(0, 1);
      }
   }

   FileName = replace_filename_format_specifiers(FileName);
   if (FileName.empty())
      return std::ofstream();

   std::ofstream Out(FileName, std::ios_base::out |
                              (Append ? std::ios_base::app : std::ios_base::trunc));
   if (BlankLine && Out.good())
      Out << '\n';
   return Out;
}

#endif
