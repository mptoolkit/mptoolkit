// -*- C++ -*- $Id$

/*
  su2.cpp

  defines quantum number SU(2).

  The SU(2) name is registered statically, it is only necessary to link against this file, and
  include the header su2.h in the file that makes use of the symmetry (which ensures that
  there are no initialization order problems).
*/

#include "su2.h"
#include "symmetrybase.h"
#include "common/convertstring.h"

#include <iostream>

namespace QuantumNumbers
{

//
// SU2
//

SU2::SU2(std::string const& s)
  : j(convert_string<half_int>(s.begin(), s.end()))
{
}

void SU2::Register()
{
   RegisterStaticSymmetry<SU2>();
}

//
// Sz
//

Sz::Sz(std::string const& s)
  : m(convert_string<half_int>(s.begin(), s.end()))
{
}

} // namespace QuantumNumbers
