// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// pstream/packformat.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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
  Pack and unpack of intrinsic types to/from a buffer according to the specified format.

  Created 2004-01-17 Ian McCulloch
*/

#if !defined(PACKFORMAT_H_DHFHFJ43YR43YY4YR4Y)
#define PACKFORMAT_H_DHFHFJ43YR43YY4YR4Y

#include "binaryformats.h"

template <typename Format, typename Buffer, typename T>
inline
void pack_format(Buffer& B, T x);

template <typename Format, typename Buffer, typename T>
inline
void unpack_format(Buffer& B, T& x);

#include "packformat.cc"

#endif
