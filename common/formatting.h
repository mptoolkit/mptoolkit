// -*- C++ -*- $Id$
//
// formatting.h
//
// common formatting functions

#if !defined(MPTOOLKIT_COMMON_FORMATTING_H)
#define MPTOOLKIT_COMMON_FORMATTING_H

#include <string>
#include <complex>
#include <sstream>

// utility function - this should be somewhere else
inline
std::string
format_complex(std::complex<double> const& c)
{
   std::ostringstream Out;
   Out.precision(16);
   Out << c.real();
   if (c.imag() > 0)
      Out << '+' << c.imag() << 'i';
   else if (c.imag() < 0)
      Out << c.imag() << 'i';
   Out.flush();
   return Out.str();
}

#endif
