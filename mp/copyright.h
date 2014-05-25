// -*- C++ -*- $Id$

#if !defined(COPYRIGHT_H_HCKJEWHFIUEWHIUGFIL)
#define COPYRIGHT_H_HCKJEWHFIUEWHIUGFIL

#include <iostream>
#include <iomanip>
#include "config.h"
#include <boost/version.hpp>

#define AS_STRING(X) AS_STRING2(X)
#define AS_STRING2(X) #X

#if defined(SVN_VERSION)
#define VERSION_ PACKAGE_VERSION " (subversion tree rev " AS_STRING(SVN_VERSION) ")"
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
      "Using Boost version " << (BOOST_VERSION / 100000) 
       << "." << (BOOST_VERSION / 100 % 1000)
       << "." << (BOOST_VERSION % 100) << "\n"
      "Copyright (c) Ian McCulloch 1999-2014 All Rights Reserved\n"
      "For license conditions email " PACKAGE_BUGREPORT "\n"
      ;
}

#undef AS_STRING
#undef AS_STRING2

#endif
