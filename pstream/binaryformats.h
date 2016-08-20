// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// pstream/binaryformats.h
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
  Defines the representations for each possible binary format,
  including the representation and endianness.

  Created 2002-04-02 Ian McCulloch

  Synopsis:

  For each format, a struct is defined

  struct <format>
  {
     static const EndianType Endianness = <implementation defined>;

     template <typename Builtin>
     struct TypeTraits
     {
        typedef CurrentFormat::<binaryformat> type;
     };
  };

  This defines the mapping of each builtin type onto a binary format,
  for the given representation.  eg, for Linux_i386,
  'long' is 32 bits so Linux_i386::TypeTraits<long>::type
  is a typedef to CurrentFormat::int32.  This means that if the host platform
  is say, alpha/tru64, and we are doing I/O to a x86/linux
  formatted stream, then a long is 32-bits in the stream but
  64-bits natively.  Thus extracting a long from the stream
  widens it from 32 bits to 64 bits, and inserting a long to the
  stream narrows it from 64 bits to 32.
*/

#if !defined(BINARYFORMATS_H_DHF6R7846743YR74TRWQGYUT2)
#define BINARYFORMATS_H_DHF6R7846743YR74TRWQGYUT2

//
// include all possible format specification files here
//

#include "binaryformat_LE_L32.h"
#include "binaryformat_LE_L64.h"
#include "binaryformat_XDR.h"

#endif
