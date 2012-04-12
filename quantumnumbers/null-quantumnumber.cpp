// -*- C++ -*- $Id$

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

} // namespace QuantumNumbers
