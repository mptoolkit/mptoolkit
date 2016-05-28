// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// pstream/endian.h
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
  Calculate the endianness of an integral type

  Created 2002-04-02 Ian McCulloch
*/

#if !defined(ENDIAN_H_ERUYEUIRY4657436ETR34656ET236T)
#define ENDIAN_H_ERUYEUIRY4657436ETR34656ET236T

#include <limits>

//
// enum EndianType
//
// Represents the possible endian status of a type.
//
// NeutralEndian: no endianness (ie. size == 1)
// BigEndian:     most significant byte is at the lowest address
// LittleEndian:  least significant byte is at the lowest address
// MixedEndian:   the byte ordering is neither big nor little.
//

enum EndianType { NeutralEndian, BigEndian, LittleEndian, MixedEndian };

char const* const EndianTypeStr[] = { "NeutralEndian", "BigEndian", "LittleEndian", "MixedEndian" };

template <class T>
EndianType GetEndianType()
{
   // construct a value of type T where the bytes, from most significant to least significant,
   // are 1, 2, ..., sizeof(T)
   int Size = sizeof(T);
   T x = 0;
   for (int i = 0; i < Size; ++i)
   {
      x <<= std::numeric_limits<unsigned char>::digits;
      x += (i+1);
   }

   // Now through an alias as unsigned char*, compare the bytes in successive memory
   // locations with what we would expect for both big endian and little endian
   // and see which (if any) match up.
   unsigned char const* CharAlias = reinterpret_cast<unsigned char const*>(&x);
   bool isBigEndian = true, isLittleEndian = true;
   for (int i = 0; i < Size; ++i)
   {
      if (CharAlias[i] != (i+1)) isBigEndian = false;
      if (CharAlias[i] != (Size-i)) isLittleEndian = false;
   }
   if (isBigEndian && isLittleEndian) return NeutralEndian;
   if (isBigEndian) return BigEndian;
   if (isLittleEndian) return LittleEndian;
   return MixedEndian;
}

#endif
