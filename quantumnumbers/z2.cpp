// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// quantumnumbers/z2.cpp
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

#include "z2.h"
#include "symmetrybase.h"

namespace QuantumNumbers
{

Z2::Z2(std::string const& s)
  : x(convert_string<int>(s.begin(), s.end()))
{
}

std::string
Z2::ToString() const
{
   return ConvertToString(x);
}

void Z2::Register()
{
   RegisterStaticSymmetry<Z2>();
}

} // namespace QuantumNumbers
