// -*- C++ -*- $Id$

#include <algorithm>

// despite how ugly this looks, it should be prety fast.

struct IsDotOrSlash
{
   inline bool operator()(char c) const { return c == '.' || c == '/'; }
};

template <class FwdIter>
half_int 
convert_string_partial<half_int, FwdIter>::apply(FwdIter start, FwdIter end)
{
   // get the integral part of the string.  This is delimited by either the end of the string, or '.' or '/'
   FwdIter end_int = std::find_if(start, end, IsDotOrSlash());

   // if this is the complete string, then its an integer
   if (end_int == end) return convert_string<int>(start, end);

   // If we encountered a '/' then the number should be of the form N/2
   if (*end_int == '/')
   {
      int n = convert_string<int>(start, end_int);
      ++end_int;
      if (end_int == end) throw string_conversion_not_enough_chars();
      if (*end_int != '2') throw string_conversion_invalid_char();
      if (++end_int != end) throw string_conversion_too_many_chars();
      return half_int(n, half_int::twice_tag());
   }

   // else *end_int == '.'
   int sign;
   if (*start == '-')
   {
      ++start;
      sign = -1;
   }
   else sign = 1;

   int n = convert_string<int>(start, end_int) * 2;

   ++end_int;
   if (end_int == end)
   {
      return half_int(n * sign, half_int::twice_tag());
   }
   else if (*end_int == '0')
   {
      // check for no extraneous characters
      if (++end_int != end) throw string_conversion_too_many_chars();
      return half_int(n * sign, half_int::twice_tag());
   }
   else if (*end_int != '5') throw string_conversion_invalid_char();
   else if (++end_int != end) throw string_conversion_too_many_chars();

   // fall-through case has last digit = '5'
   return half_int((n+1) * sign, half_int::twice_tag());  // n+1 here for the '.5'
}
