// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// quantumnumbers/su2.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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
