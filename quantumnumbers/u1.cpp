// -*- C++ -*- $Id$

#include "u1.h"
#include "symmetrybase.h"

namespace QuantumNumbers
{

U1::U1(std::string const& s)
  : x(convert_string<half_int>(s.begin(), s.end()))
{
}

std::string 
U1::ToString() const
{
   return ConvertToString(x);
}

void U1::Register()
{
   RegisterStaticSymmetry<U1>();
}

#if 0
std::string 
U1Projection::ToString() const
{
   return std::string();
}
#endif

} // namespace QuantumNumbers
