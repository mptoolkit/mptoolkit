// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// quantumnumbers/null-quantumnumber.cpp
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
  null-quantumnumber.cpp
*/

#include "null-quantumnumber.h"
#include "symmetrybase.h"

namespace QuantumNumbers
{

NullQN::NullQN(std::string const& s)
{
}

std::string
NullQN::ToString() const
{
   return "0";
}

void NullQN::Register()
{
   RegisterStaticSymmetry<NullQN>();
}

SymmetryList NullSymmetryList("Null:Null");

QuantumNumber Null(NullSymmetryList);

} // namespace QuantumNumbers
