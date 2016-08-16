// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// attic/metautility.h
//
// Copyright (C) 2002-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
