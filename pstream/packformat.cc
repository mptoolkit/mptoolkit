// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// pstream/packformat.cc
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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

#include "packendian.h"
#include "formatcast.h"

namespace Private
{

template <typename FormatType, EndianType FormatEndian>
struct PackFormatHelper
{
   template <typename Buffer, typename Type>
   inline
   static void pack(Buffer& B, Type const& x)
   {
      pack_endian<FormatEndian>(B, cast_builtin_to_format<FormatType>(x));
   }

   template <typename Buffer, typename Type>
   inline
   static void unpack(Buffer& B, Type& x)
   {
      FormatType t;
      unpack_endian<FormatEndian>(B, t);
      x = cast_format_to_builtin<Type>(t);
   }

};

} // namespace Private

template <typename Format, typename Buffer, typename T>
inline
void pack_format(Buffer& B, T x)
{
   Private::PackFormatHelper<typename Format::template TypeTraits<T>::type,
      Format::Endianness>::pack(B, x);
}

template <typename Format, typename Buffer, typename T>
inline
void unpack_format(Buffer& B, T& x)
{
   Private::PackFormatHelper<typename Format::template TypeTraits<T>::type,
      Format::Endianness>::unpack(B, x);
}
