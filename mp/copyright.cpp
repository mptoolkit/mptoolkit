// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/copyright.cpp
//
// Copyright (C) 2002-2024 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2022 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#include "copyright.h"

#if !defined(HAVE_CONFIG_H)
#error "need config.h to proceed!"
#endif
#include "config.h"

#include "common/blas_vendor.h"
#include <iostream>
#include <iomanip>
#include <boost/version.hpp>
#include <boost/algorithm/string.hpp>
#include <time.h>

#include <complex>
#include <sstream>

#define AS_STRING(X) AS_STRING2(X)
#define AS_STRING2(X) #X

#if defined(SVN_VERSION)
#define VERSION_ PACKAGE_VERSION " (subversion tree rev " AS_STRING(SVN_VERSION) ")"
#elif defined(GIT_VERSION)
#define VERSION_ PACKAGE_VERSION " (git version " AS_STRING(GIT_VERSION) ")"
#else
#define VERSION_ PACKAGE_VERSION
#endif

#if !defined(NDEBUG)
#define VERSION VERSION_ " (DEBUG)"
#else
#define VERSION VERSION_
#endif

#if defined(HAVE_OPENMP)
#define OPENMP_INFO "Using OpenMP with compiler flags " CONFIG_OPENMP_CXXFLAGS "\n"
#else
#define OPENMP_INFO
#endif

#define STRINGIZE(x) #x

void print_copyright(std::ostream& out)
{
   out << "Matrix Product Toolkit version " VERSION "\n"
      "Copyright (C) Matrix Product Toolkit Developers 1999-2024\n"
      "Compiled on " __DATE__ " at " __TIME__ "\n"
      "Configured using compiler " CONFIG_COMPILER_VENDOR " " CONFIG_COMPILER_VERSION "\n"
      "Compiler flags: " CONFIG_CXXFLAGS "\n"
      OPENMP_INFO
      // "Using Boost version " << (BOOST_VERSION / 100000) << "." << (BOOST_VERSION / 100 % 1000) << "." << (BOOST_VERSION % 100) << "\n"
      "BLAS vendor: " << BLAS_Vendor() << "\n"
      "BLAS version: " << BLAS_Version() << "\n"
      "Post bugs at: " PACKAGE_BUGREPORT "\n"
      "This program comes with ABSOLUTELY NO WARRANTY; for details run 'mp-info --warranty'.\n"
      "This is free software, and you are welcome to redistribute it under certain conditions;\n"
      "run 'mp-info --copying' for details.\n"
      "Research publications making use of this software should include appropriate citations\n"
      "and/or acknowledgements; run 'mp-info --citations' for details.\n"
      "Website: https://mptoolkit.qusim.net/\n"
      ;

   if (BLAS_Vendor() == "MKL")
   {
      std::string v = BLAS_Version();
      if (v.find("2020.0.4") != std::string::npos)
      {
         out << std::flush;
         std::cerr <<  "\n***********************************************************************************\n"
                         "*  WARNING   WARNING   WARNING   WARNING   WARNING   WARNING   WARNING   WARNING  *\n"
                         "*  The MKL version is 2020.0.4, which contains bugs in the SVD. Do not use!       *\n"
                         "*  WARNING   WARNING   WARNING   WARNING   WARNING   WARNING   WARNING   WARNING  *\n"
                         "***********************************************************************************\n\n";
         std::exit(EXIT_FAILURE);
      }
   }
}

void print_warranty(std::ostream& out)
{
   out << "\
  THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY\n\
APPLICABLE LAW.  EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT\n\
HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM \"AS IS\" WITHOUT WARRANTY\n\
OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO,\n\
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR\n\
PURPOSE.  THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM\n\
IS WITH YOU.  SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF\n\
ALL NECESSARY SERVICING, REPAIR OR CORRECTION.\n\
";
}

void print_copying(std::ostream& out)
{
   out << "This program is distributed according to the terms of the GPLv3.\n"
      "See the file COPYING in the source directory for full license conditions.\n"
      "Research publications making use of this software, including modified versions,\n"
      "should include appropriate citations and/or acknowledgements, as described in the\n"
      "file CITATIONS in the source directory, or run 'mp-info --citations'.\n";
}

void print_citations(std::ostream& out)
{
   out << "See https://github.com/mptoolkit\n";
}

#undef AS_STRING
#undef AS_STRING2

std::string Wikify(std::string const& x, bool itool)
{
   std::string Result;
   bool Capital = true;
   bool Hyphen = false;
   for (char c : x)
   {
      if (c == '-')
      {
         Hyphen = true;
         Capital = true;
      }
      else
      {
         if (Capital)
         {
            Result += char(toupper(c));
            if (!Hyphen || (!itool) || c != 'i')
               Capital = false;
            Hyphen = false;
         }
         else
            Result += char(tolower(c));
      }
   }
   return Result;
}

void print_copyright(std::ostream& out, std::string const& Category, std::string const& Name)
{
   print_copyright(out);
   std::string URL = Wikify(Category) + "/" + Wikify(Name, Name != "mp-info");
   out << "Documentation: "
       << " https://mptoolkit.qusim.net/" << URL << "\n";
}

// The basename() function is useful in help messages for printing the program name
std::string basename(std::string const& FName)
{
   return std::string(boost::find_last(FName, "/").begin(), FName.end());
}

// Escape an argument for bash (eg, to print the command-line arguments)
std::string
EscapeArgument(std::string const& s)
{
   if (s.find_first_of(" |*#^&;<>\n\t(){}[]$\\`'") != std::string::npos)
   {
      std::string Result = s;
      // escape some special characters explicitly
      boost::algorithm::replace_all(Result, "\\", "\\\\");
      boost::algorithm::replace_all(Result, "\"", "\\\"");
      boost::algorithm::replace_all(Result, "$", "\\$");
      return '"' + Result + '"';
   }
   else
      return s;
}

std::string EscapeCommandline(int argc, char** argv)
{
   std::string s = EscapeArgument(argv[0]);
   for (int i = 1; i < argc; ++i)
      s = s + ' ' + EscapeArgument(argv[i]);
   return s;
}

// Print a useful preamble, consisting of the
// program name and arguments, and the date.
void print_preamble(std::ostream& out, int argc, char** argv)
{
   out << "#" << EscapeArgument(argv[0]);
   for (int i = 1; i < argc; ++i)
      out << ' ' << EscapeArgument(argv[i]);
   time_t now = time(NULL);
   char s[200];
   int const max = 200;
   strftime(s, max, "%a, %d %b %Y %T %z", localtime(&now));
   out << "\n#Date: " << s;
   out << std::endl;
}
