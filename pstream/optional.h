/* -*- C++ -*- $Id$

  Seralization of boost::optional
*/

#if !defined(OPTIONAL_H_JSDCHWUI4357784Y7WEHOLWEHO)
#define OPTIONAL_H_JSDCHWUI4357784Y7WEHOLWEHO

#include "pstream.h"
#include <boost/optional.hpp>

namespace PStream
{

template <int format, typename T>
opstreambuf<format>& operator<<(opstreambuf<format>& out, boost::optional<T> const& x)
{
   bool b = x;
   out << b;
   if (b)
      out << (*x);
   return out;
}

template <int format, typename T>
ipstreambuf<format>& operator>>(ipstreambuf<format>& in, boost::optional<T>& x)
{
   bool b; in >> b;
   if (b)
   {
      T n;
      in >> n;
      x = n;
   }
   else
      x = boost::optional<T>();
   return in;
}

} // namespace PStream

#endif
