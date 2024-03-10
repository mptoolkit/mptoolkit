// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// pheap/polycast.h
//
// Copyright (C) 2002-2016 Ian McCulloch <ian@qusim.net>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

/*
  Type-safe, template-safe and extensible version of dynamic_cast.

  Created 2002-12-07 Ian McCulloch

  To use: exactly the same as dynamic_cast, ie
  T* t = poly_cast<U*>(u);

  The intention is that poly_cast also works for smart pointers:
  smart_ptr<T> t = poly_cast<smart_ptr<t> >(u);

  To implement poly_cast for a new type, make a (partial) specialization of
  the structure poly_cast_helper<T, U>, with a static member apply(u) which
  takes a single parameter of type U and returns a value of type T.
*/

#if !defined(POLYCAST_H_HFUIOWHF4389YHEWPHJWFE8HP)
#define POLYCAST_H_HFUIOWHF4389YHEWPHJWFE8HP

template <class T, class U>
struct poly_cast_helper;

template <class T, class U>
inline
T poly_cast(U const& val)
{
   return poly_cast_helper<T,U>::apply(val);
}

// specialization for pointer types
template <class T, class U>
struct poly_cast_helper<T*, U*>
{
   static T* apply(U* u) { return dynamic_cast<T*>(u); }
};

#endif
