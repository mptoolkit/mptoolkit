// -*- C++ -*- $Id$

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
