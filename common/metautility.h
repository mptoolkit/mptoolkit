// -*- C++ -*- $Id$

/*
  metautility.h

  metaprogramming utility functions.

  Created 2002-04-09 Ian McCulloch
*/

#if !defined(METAUTILITY_H_EUIYR7834674T7R83478O47843783O4)
#define METAUTILITY_H_EUIYR7834674T7R83478O47843783O4

// traits type to determine if two types are equal

template <class T, class U>
struct compare_types
{
   static bool const equal = false;
   static bool const not_equal = true;
};

template <class T>
struct compare_types<T, T>
{
   static bool const equal = true;
   static bool const not_equal = false;
};

// compile-time version of operator ?:

template <bool Selector, typename T1, typename T2>
struct select_type
{
   typedef T1 type;
};

template <typename T1, typename T2>
struct select_type<false, T1, T2>
{
   typedef T2 type;
};

#endif
