// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// attic/ctassert.h
//
// Copyright (C) 2016 Ian McCulloch <ian@qusim.net>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Reseach publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

/*
  ctassert.h

  Compile-time assert templates.

  Created 2000-02-28 Ian McCulloch

  Last modified 2000-05-05, renamed compile_time_assert to ct_assert

  Provides template <bool> ct_assert.  This
  class produces a compile time error if it is instantiated
  with a value of false.

  template <class T, class U> compare_types is a traits
  class that provides static const bool equal, which
  is true if T and U are the same type, and false otherwise.

  template <class T, class U> assert_equal_types
  produces a compile time error if T and U are not equal.

  usage:

  template <class Container1, class Container2>
  void foo(Container1& c1, Container2& c2)
  {

     // produce a compile-time error if c1 and c2 contain different types
     assert_equal_types<typename c1::value_type, typename c2::value_type>();

     // run different code depending on whether c1 and c2 contain different types
     // (gotta remember both if branches need to be compilable though)
     if (compare_types<typename c1::value_type, typename c2::value_type>::equal)
     {
         // code to run if c1::value_type == c2::value_type
     }
     else
     {
         //code to run if c1::value_type != c2::value_type
     }

     // assert some condition that can be evaluated at compile time
     ct_assert<sizeof(Container1::pointer) == sizeof(void*)>();
   }

   typedef usage:

   typedef ct_assert<expression> dummy_type;

*/

#if !defined(CTASSERT_H_FD7853D7ENH347Y23RY7FHUI34TY78)
#define CTASSERT_H_FD7853D7ENH347Y23RY7FHUI34TY78

#include "metautility.h"

template <bool b>
struct ct_assert
{
   int assert_failed[int(b)];
};

template <>
struct ct_assert<true> { };

template <class T, class U>
struct assert_equal_types : public ct_assert<compare_types<T, U>::equal>
{
};

#endif
