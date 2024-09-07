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

#if !defined(MPTOOLKIT_COMMON_ENVIRONMENT_H)
#define MPTOOLKIT_COMMON_ENVIRONMENT_H

#include <string>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <ctime>
#include <cstdlib>
#include <map>
#include <regex>
#include <iomanip>
#include <fstream>

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
      char const* Str = std::getenv(str.c_str());
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

inline
char const*
getenv_or_default(std::string const& str, char const* Default)
{
   char const* Str = std::getenv(str.c_str());
   return Str ? Str : Default;
}

// Format the current date/time according to the given format specifier, using strftime
// Could be replaced by Chrono in C++17
inline
std::string format_date_time(const std::string& format)
{
   std::time_t t = std::time(nullptr);
   std::tm tm = *std::localtime(&t);
   char buffer[100];
   std::strftime(buffer, sizeof(buffer), format.c_str(), &tm);
   return std::string(buffer);
}

// Find an unused filename of the form base_name N extension, where N is an integer, zero padded to minimum width
inline
std::string find_next_unused_filename(const std::string& base_name, const std::string& extension, int width = 0)
{
   int number = 1;
   std::string filename;
   std::ostringstream oss;

   // Loop until an unused filename is found
   do {
      oss.str(""); // Clear the stream
      oss.clear(); // Clear the error state
      oss << std::setw(width) << std::setfill('0') << number++;
      filename = base_name + oss.str() + extension;
   } while (std::ifstream(filename)); // Check if file exists
   return filename;
}

// Given a filename that might contain some format specifiers, substitute the values.
inline
std::string replace_format_specifiers(std::string filename)
{
   // Format mapping.  Add more here as required
   std::map<std::string, std::string> format_map = {
      {"%D", "%Y-%m-%d-%H:%M:%S"},    // ISO date-time
      {"%T", "%H:%M:%S"},             // Time only
      {"%F", "%F"},                   // year-month-date,  equivalent to %Y-%m-%d
      {"%C", "%Y%m%d_%H%M%S"},        // compact timestamp
      {"%z", "%z"},                   // timezone
      {"%Z", "%Z"}                    // timezone name
   };

   // Replace date/time format specifiers
   for (auto const& p : format_map)
   {
      std::size_t pos;
      while ((pos = filename.find(p.first)) != std::string::npos)
      {
         filename.replace(pos, p.first.length(), format_date_time(p.second));
      }
   }

   // Handle %n with optional width specifier using regex
   std::regex re(R"((.*?)(%(\d*)n)(.*))");
   std::smatch match;
   if (std::regex_match(filename, match, re))
   {
      std::string base_name = match[1].str(); // Part before %n
      // std::string full_n_spec = match[2].str(); // Full %n specifier
      std::string width_str = match[3].str(); // Digits captured
      std::string extension = match[4].str(); // Part after %n

      int width = width_str.empty() ? 0 : std::stoi(width_str); // Default width is 3 if unspecified

      // Generate a new filename with the next available number
      filename = find_next_unused_filename(base_name, extension, width);
   }

   return filename; // No %n found
}

// Open a file using the given filename template
inline
std::ofstream open_bench_file(std::string filename)
{
   if (filename.empty())
      return std::ofstream();  // Return an empty ofstream if no filename provided

   // Check for append mode: '+' at the start of the filename
   bool append_mode = (filename[0] == '+');
   bool blank = false;
   if (append_mode)
   {
      filename = filename.substr(1); // Remove the '+' sign in append mode
      if (!filename.empty() && filename[0] == '+')
      {
         blank = true;
         filename = filename.substr(1); // Remove the additional '+' sign for blank line
      }
   }

   // Replace format specifiers
   filename = replace_format_specifiers(filename);

   // Open the file
   std::ofstream out(filename, std::ios_base::out | (append_mode ? std::ios_base::app : std::ios_base::trunc));
   if (blank)
      out << '\n';
   return out;
}

#endif
