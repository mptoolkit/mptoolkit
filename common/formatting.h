// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// common/formatting.h
//
// Copyright (C) 2015-2020 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2021 Jesse Osborne <j.osborne@uqconnect.edu.au>
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
// formatting.h
//
// common formatting functions

#if !defined(MPTOOLKIT_COMMON_FORMATTING_H)
#define MPTOOLKIT_COMMON_FORMATTING_H

#include <string>
#include <complex>
#include <sstream>
#include <iomanip>

namespace formatting
{

inline
std::string
format_complex(std::complex<double> const& c)
{
   std::ostringstream Out;
   Out.precision(16);
   Out << c.real();
   if (c.imag() > 0)
      Out << "+" << c.imag() << 'i';
   else if (c.imag() < 0)
      Out << c.imag() << 'i';
   Out.flush();
   return Out.str();
}

inline
int digits(double x)
{
   using std::abs;
   using std::floor;
   using std::pow;
   double const eps = std::numeric_limits<double>::epsilon() * 8;
   int d = 0;
   while (abs(x - floor(x * pow(10,d)) / pow(10,d)) > eps)
   {
      ++d;
   }
   return d;
}

inline
int digits(std::complex<double> x)
{
   using std::real;
   using std::imag;
   using std::max;
   return max(digits(real(x)), digits(imag(x)));
}

inline
std::string
format_digits(double x, int Digits)
{
   std::ostringstream str;
   str << std::fixed << std::setprecision(Digits) << std::setfill('0') << x;
   return str.str();
}

inline
std::string
format_digits(std::complex<double> c, int Digits)
{
   std::string Result = format_digits(c.real(), Digits);
   if (c.imag() > 0)
      Result += "+" + format_digits(c.imag(), Digits) + 'i';
   else if (c.imag() < 0)
      Result += format_digits(c.imag(), Digits) + 'i';
   return Result;
}

} // namespace formatting

#endif
