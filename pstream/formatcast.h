// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// pstream/formatcast.h
//
// Copyright (C) 2003-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
  template function cast_builtin_to_format<T>(src) to cast a builtin type to a format type,
  and template function cast_format_to_builtin<T>(src) to cast a format type to a builtin.

  Created 2003-03-04 Ian McCulloch

  This is basically a wrapper around boost::numeric_cast, with some optimizations
  for trivial conversions that boost does not yet do (and look hard for compilers to optimize).

  Also, it is guaranteed that all format types are supported as source or destination, which
  isn't necessarily true of boost::numeric_cast.
*/

#if !defined(FORMATCAST_H_SADKJEHRU387RY43875YFHO9HE8F4)
#define FORMATCAST_H_SADKJEHRU387RY43875YFHO9HE8F4

#include <boost/cast.hpp>

namespace Private
{

template <typename DestType, bool isTrivial>
  struct FormatCastHelper  // default version handles isTrivial == false
{
   template <typename SrcType>
   static DestType apply(SrcType f) { return boost::numeric_cast<DestType>(f); }
};

template <typename DestType>
struct FormatCastHelper<DestType, true>
{
   template <typename SrcType>
   static DestType apply(SrcType f) { return static_cast<DestType>(f); }
};

} // namespace Private

template <typename FormatType, typename BuiltinType>
inline
FormatType cast_builtin_to_format(BuiltinType f)
{
   return Private::FormatCastHelper<FormatType,
     CurrentFormat::IsTrivialConversion<BuiltinType, FormatType>::value>::apply(f);
}

template <typename BuiltinType, typename FormatType>
inline
BuiltinType cast_format_to_builtin(FormatType f)
{
   return Private::FormatCastHelper<BuiltinType,
     CurrentFormat::IsTrivialConversion<BuiltinType, FormatType>::value>::apply(f);
}

#endif
