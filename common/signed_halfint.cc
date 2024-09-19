// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// common/signed_halfint.cc
//
// Copyright (C) 2016 Ian McCulloch <ian@qusim.net>
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

// inlines

inline
bool operator==(signed_half_int h1, signed_half_int h2)
{
   return h1.twice() == h2.twice();
}

inline
bool operator!=(signed_half_int h1, signed_half_int h2)
{
   return h1.twice() != h2.twice();
}

inline
bool operator<(signed_half_int h1, signed_half_int h2)
{
   return h1.twice() < h2.twice();
}

inline
bool operator>(signed_half_int h1, signed_half_int h2)
{
   return h1.twice() > h2.twice();
}

inline
bool operator<=(signed_half_int h1, signed_half_int h2)
{
   return h1.twice() <= h2.twice();
}

inline
bool operator>=(signed_half_int h1, signed_half_int h2)
{
   return h1.twice() >= h2.twice();
}

inline
half_int abs(signed_half_int h)
{
   return half_int(h.twice_magnitude(), signed_half_int::twice_tag());
}

inline
int sign(signed_half_int h)
{
   return h.sign();
}

#if defined(USE_PSTREAM)
inline
PStream::opstream& operator<<(PStream::opstream& out, const signed_half_int& H)
{
   return out << H.N2;
}

inline
PStream::ipstream& operator>>(PStream::ipstream& in, signed_half_int& H)
{
   return in >> H.N2;
}
#endif // #defined (USE_PSTREAM)
