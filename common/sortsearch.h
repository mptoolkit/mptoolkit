// -*- C++ -*- $Id$

/*
  sortsearch.h

  A collection of simple sorting and searching template functions.

  Created 2001-08-17 from implementations scattered across a few places

  The implementations of sort3, sort4 are not very good.
*/

#if !defined(SORTSEARCH_H_FH34789Y348UYJFIO45897678ERYUI34HIOWEJ7Y)
#define SORTSEARCH_H_FH34789Y348UYJFIO45897678ERYUI34HIOWEJ7Y

#include <algorithm>

// generalizations of min() and max() to larger numbers of arguments

template <class T>
inline
T max3(T const& a, T const& b, T const& c)
{
   using std::max;
   return max(max(a,b), c);
}

template <class T>
inline
T min3(T const& a, T const& b, T const& c)
{
   using std::min;
   return min(min(a,b), c);
}

template <class T>
inline
T max4(T const& a, T const& b, T const& c, T const& d)
{
   using std::max;
   return max(max3(a,b,c),d);
}

template <class T>
inline
T min4(T const& a, T const& b, T const& c, T const& d)
{
   using std::min;
   return min(min3(a,b,c),d);
}

template <class T>
inline
T max5(T const& a, T const& b, T const& c, T const& d, T const& e)
{
   return max3(max3(a,b,c),d,e);
}

template <class T>
inline
T min5(T const& a, T const& b, T const& c, T const& d, T const& e)
{
   return min3(min3(a,b,c),d,e);
}

// sorting a fixed number of arguments into increasing order using operator<

template <class T>
inline 
bool sort2(T& a, T& b)
{
   using std::swap;
   if (b < a) 
   {
      swap(a,b);
      return false;
   }
   return true;
}

template <class T>
inline
void sort3(T& a, T& b, T& c)
{
   using std::swap;
   if (b < a)
   {
      if (c < b)
      {
	 // order is c < b < a
	 swap(a,c);
      }
      else
      {
	 // b <= c
	 swap(a,b);
	 if (c < b) swap(b,c);
      }
   }
   else
   {
      // a <= b
      if (c < b)
      {
	 swap(b,c);
	 if (b < a) swap(a,b);
      }
   }
}

template <class T>
inline 
void sort4(T& a, T& b, T& c, T& d)
{
   // this isn't as efficient as it could be
   sort3(a,b,c);
   sort2(c,d);
   sort2(b,c);
   sort2(a,b);
}

#endif


