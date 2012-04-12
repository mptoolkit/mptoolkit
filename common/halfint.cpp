// -*- C++ -*- $Id$

#include "halfint.h"
#include <iomanip>

std::ostream& operator<<(std::ostream& out, const half_int& H)
{
   out << (H.twice() / 2.0);
   return out;
}

std::istream& operator>>(std::istream& in, half_int& H)
{
#if 0
   // this is a (possibly) faster (except for the string allocation), but less robust version
   in >> skipws;
   char c;
   std::string s;
   in >> c;
   s = c;
   TRACE(c);
   while ((in >> c) && ((c >= '0' && c <= '9') || c == '.' || c == '/'))
   {
      s += c;
   }
   in.putback(c);
   TRACE(s);
   H = convert_string<half_int>(s.begin(), s.end());
   return in;
#else
   double d;
   in >> d; 
   H = half_int(d);
   return in;
#endif
}

void half_int::throw_cannot_convert()
{
   throw std::runtime_error("half_int: cannot convert odd half_int to integer!"); 
}
