// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// pstream/binaryformat_LE_L32.h
//
// Copyright (C) 2002-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

/*
  Configuration file for binary format LE_L32

  Created 2002-04-02 Ian McCulloch

  The format LE_L32 is little-endian, with
  bool      -> 8 bits
  char      -> 8 bits
  short     -> 16 bits
  int       -> 32 bits
  long      -> 32 bits
  long long -> 64 bits
  float     -> 32 bit IEEE
  double    -> 64 bit IEEE

  This is the native format for Linux on x86 architecture.
*/

#if !defined(BINARYFORMAT_LE_L32_H_KJSHFIUEWRI43H5Y7587Y87OHYY187HIUEWL)
#define BINARYFORMAT_LE_L32_H_KJSHFIUEWRI43H5Y7587Y87OHYY187HIUEWL

#include "currentformat.h"

namespace PStream
{

namespace BinaryFormats
{

struct LE_L32
{
   static const EndianType Endianness = LittleEndian;

   typedef CurrentFormat::uint32 size_type;
   typedef CurrentFormat::int32  difference_type;

   template <class T>
   struct TypeTraits;
};

template <>
struct LE_L32::TypeTraits<bool>
{
   typedef CurrentFormat::char8 type;
};

template <>
struct LE_L32::TypeTraits<char>
{
   typedef CurrentFormat::char8 type;
};

template <>
struct LE_L32::TypeTraits<signed char>
{
   typedef CurrentFormat::char8 type;
};

template <>
struct LE_L32::TypeTraits<unsigned char>
{
   typedef CurrentFormat::uchar8 type;
};

template <>
struct LE_L32::TypeTraits<short>
{
   typedef CurrentFormat::int16 type;
};

template <>
struct LE_L32::TypeTraits<unsigned short>
{
   typedef CurrentFormat::uint16 type;
};

template <>
struct LE_L32::TypeTraits<int>
{
   typedef CurrentFormat::int32 type;
};

template <>
struct LE_L32::TypeTraits<unsigned int>
{
   typedef CurrentFormat::uint32 type;
};

template <>
struct LE_L32::TypeTraits<long>
{
   typedef CurrentFormat::int32 type;
};

template <>
struct LE_L32::TypeTraits<unsigned long>
{
   typedef CurrentFormat::uint32 type;
};

#if defined(USE_LONGLONG)
template <>
struct LE_L32::TypeTraits<long long>
{
   typedef CurrentFormat::int64 type;
};

template <>
struct LE_L32::TypeTraits<unsigned long long>
{
   typedef CurrentFormat::uint64 type;
};
#endif

template <>
struct LE_L32::TypeTraits<float>
{
   typedef CurrentFormat::ieee32 type;
};

template <>
struct LE_L32::TypeTraits<double>
{
   typedef CurrentFormat::ieee64 type;
};

} // namespace BinaryFormats

} // namespace PStream

#endif
