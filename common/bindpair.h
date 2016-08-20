// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/bindpair.h
//
// Copyright (C) 1999-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
  bindpair.h

  Created 10/6/99
  Ian McCulloch

  Defines the reference_pair template class, and the bind_pair() function.  The reference_pair
  class has no public constructors and should not be used directly.
  bind_pair(x, y) returns a pair-like object that acts as an l-value.  For example,

  pair<int, int> MyPair;
  int x, y;

  MyPair = make_pair(4, 5);
  bind_pair(x,y) = MyPair;
  // now x == 4 and y == 5

  We also define a couple of (unrelated) functors, extract_first, extract_second,
  that return first and second respectively of a pair.
*/

#if !defined(BINDPAIR_H_D64356UK7Y7F53FYYY755644R783478Y89UR567HUY)
#define BINDPAIR_H_D64356UK7Y7F53FYYY755644R783478Y89UR567HUY

#include <utility>
#include <functional>

template <class T1, class T2> class reference_pair;

template <class T1, class T2>
reference_pair<T1,T2> bind_pair(T1& first, T2& second);

template <class T1, class T2>
class reference_pair
{
   public:
      void operator=(const std::pair<T1, T2>& Value) { first = Value.first; second = Value.second; }

   private:
      // constructor, not for public use!
      reference_pair(T1& FirstRef, T2& SecondRef) : first(FirstRef), second(SecondRef) {}

   private:
      reference_pair(const reference_pair<T1, T2>& V) : first(V.first), second(V.second) {}

      void operator=(reference_pair const&); // not implemented

      T1& first;
      T2& second;

      // SGI compiler is broken w.r.t. friend template functions, so we have to introduce more parameters
      // template <class U1, class U2> friend reference_pair<U1,U2> bind_pair(U1& first, U2& second);

   friend reference_pair<T1,T2> bind_pair<T1, T2>(T1& first, T2& second);
};

template <class T1, class T2>
inline reference_pair<T1,T2> bind_pair(T1& first, T2& second)
{
   return reference_pair<T1,T2>(first, second);
}

template <class T1, class T2>
struct extract_first : public std::unary_function<std::pair<T1, T2>, T1>
{
   const T1& operator()(const std::pair<T1, T2>& x) const { return x.first; }
};

template <class T1, class T2>
struct extract_second : public std::unary_function<std::pair<T1, T2>, T2>
{
   const T2& operator()(const std::pair<T1, T2>& x) const { return x.second; }
};

#endif
