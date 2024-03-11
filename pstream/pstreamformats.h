// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// pstream/pstreamformats.h
//
// Copyright (C) 2002-2016 Ian McCulloch <ian@qusim.net>
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
  This defines the binary formats used by PStreams

  Created 2002-04-11 Ian McCulloch,
*/

#if !defined(PSTREAMFORMATS_H_HDJKFHUIY43Q897Y87QY9P87QAP89T3Q4Y89)
#define PSTREAMFORMATS_H_HDJKFHUIY43Q897Y87QY9P87QAP89T3Q4Y89

#include "binaryformats.h"



#include "packformat.h"
#include <iostream>
#include <typeinfo>
#include "common/trace.h"
#include "common/inttype.h"

namespace PStream
{

// The int Format values are defined in namespace PStream::format
//
// Currently implemented formats:
//
// Native      Fast but non-portable
// LE_L32      Little-endian with 32-bit 'long'.
// LE_L64      Little-endian with 64-bit 'long'.
// XDR         XDR-compatible, but 'long' is 64 bits. Big-endian.
// Linux_i386  alias for LE_L32.
// Tru64_Alpha alias for LE_L64.
// Host        alias for the format most closely resembling the host binary format.
//
// Don't use format::Native for persistent streams.  It is intended only for
// memory buffers to be read by the same process that generated them.
//

namespace format
{

int const Native = 0;

int const LE_L32 = 1;
int const Linux_i386 = 1;

int const LE_L64 = 2;
int const Tru64_Alpha = 2;

int const XDR = 3;

// Try to guess the current format
#if defined(WORDS_BIGENDIAN)
int const Host = XDR;
#elif SIZEOF_LONG == 8
int const Host = LE_L64;
#else
int const Host = LE_L32;
#endif

} // namespace format

//
// format_traits
//
// Determines the format type (from the types in namespace PStream::BinaryFormats)
// corresponding to the given int Format label.
//

template <int Format>
struct format_traits;

template <> struct format_traits<format::LE_L32> { typedef BinaryFormats::LE_L32 Format; };
template <> struct format_traits<format::LE_L64> { typedef BinaryFormats::LE_L64 Format; };
template <> struct format_traits<format::XDR>    { typedef BinaryFormats::XDR    Format; };

// For Native format, ultimately we want to do better than just emulating the host format.
template <> struct format_traits<format::Native> : public format_traits<format::Host> {};

//
// pstreambuf_traits determines size_type and difference_types from the respective format.
//

template <int Format>
struct pstreambuf_traits
{
   typedef typename format_traits<Format>::Format  format;
   typedef typename format::size_type              size_type;
   typedef typename format::difference_type        difference_type;
};

//
// Streaming of builtins
//

#define DEFINE_PBINARYBUF_BUILTIN(type, Which)                                          \
inline                                                                                  \
opstreambuf<Which>& operator<<(opstreambuf<Which>& out, type c)                         \
{                                                                                       \
   pack_format<pstreambuf_traits<Which>::format>(out, c);                               \
   return out;                                                                          \
}                                                                                       \
inline                                                                                  \
ipstreambuf<Which>& operator>>(ipstreambuf<Which>& in, type& c)                         \
{                                                                                       \
   unpack_format<pstreambuf_traits<Which>::format>(in, c);                              \
   return in;                                                                           \
}

#define DEFINE_ALL_PBINARYBUF_BUILTIN(type)     \
DEFINE_PBINARYBUF_BUILTIN(type, 0)              \
DEFINE_PBINARYBUF_BUILTIN(type, 1)              \
DEFINE_PBINARYBUF_BUILTIN(type, 2)              \
DEFINE_PBINARYBUF_BUILTIN(type, 3)

DEFINE_ALL_PBINARYBUF_BUILTIN(bool)
DEFINE_ALL_PBINARYBUF_BUILTIN(char)
DEFINE_ALL_PBINARYBUF_BUILTIN(signed char)
DEFINE_ALL_PBINARYBUF_BUILTIN(unsigned char)
DEFINE_ALL_PBINARYBUF_BUILTIN(short)
DEFINE_ALL_PBINARYBUF_BUILTIN(unsigned short)
DEFINE_ALL_PBINARYBUF_BUILTIN(int)
DEFINE_ALL_PBINARYBUF_BUILTIN(unsigned int)
DEFINE_ALL_PBINARYBUF_BUILTIN(long)
DEFINE_ALL_PBINARYBUF_BUILTIN(unsigned long)
#if defined(USE_LONGLONG)
DEFINE_ALL_PBINARYBUF_BUILTIN(long long)
DEFINE_ALL_PBINARYBUF_BUILTIN(unsigned long long)
#endif
DEFINE_ALL_PBINARYBUF_BUILTIN(float)
DEFINE_ALL_PBINARYBUF_BUILTIN(double)

#define DEFINE_PBINARYBUF_FIXEDSIZE(type, Which)                                        \
inline                                                                                  \
opstreambuf<Which>& operator<<(opstreambuf<Which>& out, inttype::StrongType<type> c)    \
{                                                                                       \
   pack_endian<pstreambuf_traits<Which>::format::Endianness>(out, c.value());           \
   return out;                                                                          \
}                                                                                       \
inline                                                                                  \
ipstreambuf<Which>& operator>>(ipstreambuf<Which>& in, inttype::StrongType<type>& c)    \
{                                                                                       \
   unpack_endian<pstreambuf_traits<Which>::format::Endianness>(in, c.value());          \
   return in;                                                                           \
}

#define DEFINE_ALL_PBINARYBUF_FIXEDSIZE(type)   \
DEFINE_PBINARYBUF_FIXEDSIZE(type, 0)            \
DEFINE_PBINARYBUF_FIXEDSIZE(type, 1)            \
DEFINE_PBINARYBUF_FIXEDSIZE(type, 2)            \
DEFINE_PBINARYBUF_FIXEDSIZE(type, 3)

DEFINE_ALL_PBINARYBUF_FIXEDSIZE(boost::int8_t)
DEFINE_ALL_PBINARYBUF_FIXEDSIZE(boost::uint8_t)
DEFINE_ALL_PBINARYBUF_FIXEDSIZE(boost::int16_t)
DEFINE_ALL_PBINARYBUF_FIXEDSIZE(boost::uint16_t)
DEFINE_ALL_PBINARYBUF_FIXEDSIZE(boost::int32_t)
DEFINE_ALL_PBINARYBUF_FIXEDSIZE(boost::uint32_t)
DEFINE_ALL_PBINARYBUF_FIXEDSIZE(boost::int64_t)
DEFINE_ALL_PBINARYBUF_FIXEDSIZE(boost::uint64_t)

#if 0

template <int Which, typename T>
inline
opstreambuf<Which>& operator<<(opstreambuf<Which>& out, inttype::StrongType<T> const& x)
{
   pack_endian<pstreambuf_traits<Which>::format::Endianness>(out, x.value());
   return out;
}

template <int Which, typename T>
inline
ipstreambuf<Which>& operator>>(ipstreambuf<Which>& in, inttype::StrongType<T> const& x)
{
   unpack_endian<pstreambuf_traits<Which>::format::Endianness>(in, x.value());
   return in;
}

#endif

} // namespace pstream

#endif
