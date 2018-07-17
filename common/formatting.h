// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/formatting.h
//
// Copyright (C) 2015-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Reseach publications making use of this software should include
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
#include <ostream>
#include <iomanip>

inline
std::string
format_complex(std::complex<double> const& c)
{
   std::ostringstream Out;
   Out.precision(16);
   Out << c.real();
   if (c.imag() > 0)
      Out << " + " << c.imag() << 'i';
   else if (c.imag() < 0)
      Out << ' ' << c.imag() << 'i';
   Out.flush();
   return Out.str();
}

template <typename T>
void
wite_format(std::ostream& out, T x)
{
   out << std::setw(6) << x;
}

inline
void
write_format(std::ostream& out, double x)
{
   out << std::setw(10) << x;
}

inline
void
write_format(std::ostream& out, std::complex<double> x)
{
   out << format_complex(x);
}

#endif
