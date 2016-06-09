// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// pstream/binaryformat_LE_L64.h
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
  Configuration file for binary format LE_L32

  Created 2002-04-02 Ian McCulloch

  The format LE_L64 is little-endian, with
  bool      -> 8 bits
  char      -> 8 bits
  short     -> 16 bits
  int       -> 32 bits
  long      -> 64 bits
  long long -> 64 bits
  float     -> 32 bit IEEE
  double    -> 64 bit IEEE

  This is the native format for Tru64 using the 64-bit ABI
*/

#if !defined(BINARYFORMAT_LE_L64_H_JKHOIUEWHROIUHROI32H98YR7843Y6YWFH)
#define BINARYFORMAT_LE_L64_H_JKHOIUEWHROIUHROI32H98YR7843Y6YWFH

#include "currentformat.h"

namespace PStream
{

namespace BinaryFormats
{

struct LE_L64
{
   static const EndianType Endianness = LittleEndian;

   typedef CurrentFormat::uint64 size_type;
   typedef CurrentFormat::int64  difference_type;
   
   template <class T>
   struct TypeTraits;
};

template <>
struct LE_L64::TypeTraits<bool>
{
   typedef CurrentFormat::char8 type;
};

template <>
struct LE_L64::TypeTraits<char>
{
   typedef CurrentFormat::char8 type;
};

template <>
struct LE_L64::TypeTraits<signed char>
{
   typedef CurrentFormat::char8 type;
};

template <>
struct LE_L64::TypeTraits<unsigned char>
{
   typedef CurrentFormat::uchar8 type;
};

template <>
struct LE_L64::TypeTraits<short>
{
   typedef CurrentFormat::int16 type;
};

template <>
struct LE_L64::TypeTraits<unsigned short>
{
   typedef CurrentFormat::uint16 type;
};

template <>
struct LE_L64::TypeTraits<int>
{
   typedef CurrentFormat::int32 type;
};

template <>
struct LE_L64::TypeTraits<unsigned int>
{
   typedef CurrentFormat::uint32 type;
};

template <>
struct LE_L64::TypeTraits<long>
{
   typedef CurrentFormat::int64 type;
};

template <>
struct LE_L64::TypeTraits<unsigned long>
{
   typedef CurrentFormat::uint64 type;
};

#if defined(USE_LONGLONG)
template <>
struct LE_L64::TypeTraits<long long>
{
   typedef CurrentFormat::int64 type;
};

template <>
struct LE_L64::TypeTraits<unsigned long long>
{
   typedef CurrentFormat::uint64 type;
};
#endif

template <>
struct LE_L64::TypeTraits<float>
{
   typedef CurrentFormat::ieee32 type;
};

template <>
struct LE_L64::TypeTraits<double>
{
   typedef CurrentFormat::ieee64 type;
};

} // namespace BinaryFormats

} // namespace PStream

#endif
