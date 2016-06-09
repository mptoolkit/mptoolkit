// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// pstream/packendian.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#if !defined(PACKENDIAN_H_DSHFIUYT437Y7LIUHYIURH8743Y7H)
#define PACKENDIAN_H_DSHFIUYT437Y7LIUHYIURH8743Y7H

#include "endian.h"
#include "currentformat.h"

template <EndianType Endian, typename Buffer, typename T>
inline
void pack_endian(Buffer& B, T const& x);

template <EndianType Endian, typename Buffer, typename T>
inline
void unpack_endian(Buffer& B, T& x);

#include "packendian.cc"

#endif
