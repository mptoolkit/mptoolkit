// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// pstream/castformat.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include <boost/cast.hpp>

template <typename DestType, typename SrcType>
inline
DestType format_cast(SrcType f)
{
   return Private::FormatCastHelper<DestType,
     CurrentFormat::IsTrivialConversion<SrcType, DestType::value>::apply(f);
}

namespace Private
{

template <typename DestType, bool isTrivial>
struct FormatCastHelper
{
   static DestType apply(SrcType f) { return boost::numeric_cast<DestType>(f); }
};

template <typename DestType>
struct FormatCastHelper<DestType, true>
{
   static DestType apply(DestType f) { return static_cast<DestType>(f); }
};

} // namespace Private
