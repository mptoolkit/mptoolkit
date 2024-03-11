// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// pstream/binaryformat_XDR.h
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
  Configuration file for platform

  Format_XDR

  Created 2002-11-12 Ian McCulloch

  In the 'traditional' XDR mapping, both int and long
  are 32 bits.  This is broken on 64-bit architectures,
  so we choose here to make long 64-bits, along with
  size_type and difference_type.
*/

#if !defined(BINARYFORMAT_XDR_H_KJFDG43Y78HDLIH87RHYLDAHY87WAL)
#define BINARYFORMAT_XDR_H_KJFDG43Y78HDLIH87RHYLDAHY87WAL

#include "currentformat.h"

namespace PStream
{

namespace BinaryFormats
{

struct XDR
{
   static const EndianType Endianness = BigEndian;

   typedef CurrentFormat::uint64 size_type;
   typedef CurrentFormat::int64  difference_type;

   template <class T>
   struct TypeTraits;
};

template <>
struct XDR::TypeTraits<bool>
{
   typedef CurrentFormat::char8 type;
};

template <>
struct XDR::TypeTraits<char>
{
   typedef CurrentFormat::char8 type;
};

template <>
struct XDR::TypeTraits<signed char>
{
   typedef CurrentFormat::char8 type;
};

template <>
struct XDR::TypeTraits<unsigned char>
{
   typedef CurrentFormat::uchar8 type;
};

template <>
struct XDR::TypeTraits<short>
{
   typedef CurrentFormat::int16 type;
};

template <>
struct XDR::TypeTraits<unsigned short>
{
   typedef CurrentFormat::uint16 type;
};

template <>
struct XDR::TypeTraits<int>
{
   typedef CurrentFormat::int32 type;
};

template <>
struct XDR::TypeTraits<unsigned int>
{
   typedef CurrentFormat::uint32 type;
};

template <>
struct XDR::TypeTraits<long>
{
   typedef CurrentFormat::int64 type;
};

template <>
struct XDR::TypeTraits<unsigned long>
{
   typedef CurrentFormat::uint64 type;
};

#if defined(USE_LONGLONG)
template <>
struct XDR::TypeTraits<long long>
{
   typedef CurrentFormat::int64 type;
};

template <>
struct XDR::TypeTraits<unsigned long long>
{
   typedef CurrentFormat::uint64 type;
};
#endif

template <>
struct XDR::TypeTraits<float>
{
   typedef CurrentFormat::ieee32 type;
};

template <>
struct XDR::TypeTraits<double>
{
   typedef CurrentFormat::ieee64 type;
};

} // namespace BinaryFormats

} // namespace PStream

#endif
