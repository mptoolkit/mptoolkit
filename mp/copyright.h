// -*- C++ -*- $Id$

#if !defined(COPYRIGHT_H_HCKJEWHFIUEWHIUGFIL)
#define COPYRIGHT_H_HCKJEWHFIUEWHIUGFIL

#include <iostream>
#include <iomanip>
#include "config.h"

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
      "Copyright (c) Ian McCulloch 1999-2010 All Rights Reserved\n"
      "For license conditions email " PACKAGE_BUGREPORT "\n"
      ;
}

#undef AS_STRING
#undef AS_STRING2

#endif
