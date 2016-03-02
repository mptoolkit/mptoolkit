// -*- C++ -*- $Id$

#if !defined(MPTOOLKIT_MP_COPYRIGHT_H)
#define MPTOOLKIT_MP_COPYRIGHT_H

#if !defined(HAVE_CONFIG_H)
#error "need config.h to proceed!"
#endif

#include "config.h"
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

inline void print_copyright(std::ostream& out)
{
   out << "Matrix Product Toolkit version " VERSION "\n"
      "Compiled on " __DATE__ " at " __TIME__ "\n"
      "Configured using compiler " CONFIG_COMPILER_VENDOR " " CONFIG_COMPILER_VERSION "\n"
      "Using Boost version " << (BOOST_VERSION / 100000) 
       << "." << (BOOST_VERSION / 100 % 1000)
       << "." << (BOOST_VERSION % 100) << "\n"
      "Copyright (c) Ian McCulloch 1999-2016 All Rights Reserved\n"
      "For license conditions email " PACKAGE_BUGREPORT "\n"
      "Documentation see http://physics.uq.edu.au/people/ianmcc/mptoolkit/\n"
      ;
}

#undef AS_STRING
#undef AS_STRING2

// The basename() function is useful in help messages for printing the program name
inline
std::string basename(std::string const& FName)
{
   return std::string(boost::find_last(FName, "/").begin(), FName.end());
}

// Escape an argument for bash (eg, to print the command-line arguments)
inline
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

inline
std::string EscapeCommandline(int argc, char** argv)
{
   std::string s = EscapeArgument(argv[0]);
   for (int i = 1; i < argc; ++i)
      s = s + ' ' + EscapeArgument(argv[i]);
   return s;
}

// Print a useful preamble, consisting of the
// program name and arguments, and the date.
inline
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

#endif
