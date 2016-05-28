// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// quantumnumbers/null-quantumnumber.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
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
