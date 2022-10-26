// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mpo/zero_mpo.h
//
// Copyright (C) 2015-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

// ZeroMPO
//
// It is convenient to have an explicit representation for an MPO that
// is empty, so that we can have a default-constructed InfiniteMPO
// that is a 'true zero' without having to worry about the distinction
// between BasicTriangularMPO and ProductMPO's.

#if !defined(MPTOOLKIT_MPO_ZERO_MPO_H)
#define MPTOOLKIT_MPO_ZERO_MPO_H

#include "pstream/pstream.h"

struct ZeroMPO {};

inline
PStream::opstream& operator<<(PStream::opstream& out, ZeroMPO const&)
{
   return out;
}

inline
PStream::ipstream& operator>>(PStream::ipstream& in, ZeroMPO&)
{
   return in;
}

inline
ZeroMPO pow(ZeroMPO, int)
{
   return ZeroMPO();
}

inline
ZeroMPO conj(ZeroMPO)
{
   return ZeroMPO();
}

inline
ZeroMPO abs(ZeroMPO)
{
   return ZeroMPO();
}

inline
ZeroMPO adjoint(ZeroMPO)
{
   return ZeroMPO();
}

inline
ZeroMPO inv_adjoint(ZeroMPO)
{
   return ZeroMPO();
}

inline
ZeroMPO dot(ZeroMPO, ZeroMPO)
{
   return ZeroMPO();
}

inline
ZeroMPO gauge_flip(ZeroMPO)
{
   return ZeroMPO();
}

#endif
