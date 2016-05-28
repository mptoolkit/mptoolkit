// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/signed_halfint.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
// Class to represent a signed half-integer, where
// 0 also has a sign.  +0 and -0 are distinct elements.
// This is used eg to label
// representations of the dihedral group.

#if !defined(MPTOOLKIT_COMMON_SIGNED_HALFINT_H)
#define MPTOOLKIT_COMMON_SIGNED_HALFINT_H

#include <iostream>
#include <stdexcept>
#include <math.h>
#include <string>
#include "halfint.h"

#if defined(USE_PSTREAM)
#include "pstream/pstream.h"
#endif

class signed_half_int
{
   public:
      signed_half_int() : N2(0) {}  // initialize to zero, rather than leave N2 unitialized.

      signed_half_int(int Sign, double D);
      signed_half_int(int Sign, int I);

      explicit signed_half_int(std::string const& s);

      // use compiler defined copy ctor, assignment and dtor
      // (user defined ctor affects performance on intel 6 compiler)

      int sign() const;

      double magnitude() const;

      int twice_magnitude() const;

      std::string to_string() const;

   private:
      int N2;

#if defined(USE_PSTREAM)
   friend PStream::opstream& operator<<(PStream::opstream& out, const signed_half_int& H);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, signed_half_int& H);
#endif

};

// non-members acting on a signed_half_int

bool operator==(signed_half_int h1, signed_half_int h2);
bool operator!=(signed_half_int h1, signed_half_int h2);

bool operator<(signed_half_int h1, signed_half_int h2);
bool operator>(signed_half_int h1, signed_half_int h2);
bool operator<=(signed_half_int h1, signed_half_int h2);
bool operator>=(signed_half_int h1, signed_half_int h2);

half_int abs(signed_half_int h);

int sign(signed_half_int h);

// conversion to string, as a fraction,
// if h is half-integral then the string is of the forn "n/2"
std::string to_string_fraction(signed_half_int h);

std::ostream& operator<<(std::ostream& out, const signed_half_int& H);
std::istream& operator>>(std::istream& in, signed_half_int& H);

#include "signed_halfint.cc"

#endif
