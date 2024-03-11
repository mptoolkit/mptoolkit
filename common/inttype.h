// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// common/inttype.h
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
  inttype.h

  (quasi-)portable fixed size integers.

  Created 2002-11-15 Ian McCulloch
*/

#if !defined(INTSIZE_H_CH879Q34Y89FYIOHF89EY8942YRJEOI)
#define INTSIZE_H_CH879Q34Y89FYIOHF89EY8942YRJEOI

#include <iostream>
#include <boost/cstdint.hpp>
#include <limits>

namespace inttype
{

//
// portable fixed-size integers.
// The _t versions are typedef's for builtin types.
// The versions without _t are strong types (ie, new types),
// basically for use in pstreams,
// so that they get the correct size irrespective of what the target
// format indicates.
//

using boost::int8_t;
using boost::uint8_t;
using boost::int16_t;
using boost::uint16_t;
using boost::int32_t;
using boost::uint32_t;
using boost::int64_t;
using boost::uint64_t;

template <typename T>
struct StrongType
{
   StrongType() {}
   StrongType(T x_) : x(x_) {}
   StrongType& operator=(StrongType const& x_) { x = x_.x; return *this; }
   StrongType& operator=(T x_) { x = x_; return *this; }

   operator T&() { return x; }
   operator T() const { return x; }

   T& value() { return x; }
   T const& value() const { return x; }

   private:
      T x;
};

typedef StrongType<int8_t>   int8;
typedef StrongType<uint8_t>  uint8;
typedef StrongType<int16_t>  int16;
typedef StrongType<uint16_t> uint16;
typedef StrongType<int32_t>  int32;
typedef StrongType<uint32_t> uint32;
typedef StrongType<int64_t>  int64;
typedef StrongType<uint64_t> uint64;

#define DECLARE_INTTYPE_STREAM_OPERATOR(type)                   \
inline std::ostream& operator<<(std::ostream& out, type x)      \
{                                                               \
   return out << type##_t(x);                                   \
}                                                               \
inline std::istream& operator>>(std::istream& in, type& x)      \
{                                                               \
   type##_t temp;                                               \
   in >> temp;                                                  \
   x = temp;                                                    \
   return in;                                                   \
}

DECLARE_INTTYPE_STREAM_OPERATOR(int8)
DECLARE_INTTYPE_STREAM_OPERATOR(int16)
DECLARE_INTTYPE_STREAM_OPERATOR(int32)
DECLARE_INTTYPE_STREAM_OPERATOR(int64)
DECLARE_INTTYPE_STREAM_OPERATOR(uint8)
DECLARE_INTTYPE_STREAM_OPERATOR(uint16)
DECLARE_INTTYPE_STREAM_OPERATOR(uint32)
DECLARE_INTTYPE_STREAM_OPERATOR(uint64)

#undef DECLARE_INTTYPE_STREAM_OPERATOR

} // namespace inttype

// Declare specializations of std::numeric_limits<> for the strong typedef's
namespace std
{

#define DEFINE_INTTYPE_LIMIT(type)                                                              \
template <> struct numeric_limits<inttype::type> : public numeric_limits<inttype::type##_t> {};

DEFINE_INTTYPE_LIMIT(int8)
DEFINE_INTTYPE_LIMIT(uint8)
DEFINE_INTTYPE_LIMIT(int16)
DEFINE_INTTYPE_LIMIT(uint16)
DEFINE_INTTYPE_LIMIT(int32)
DEFINE_INTTYPE_LIMIT(uint32)
DEFINE_INTTYPE_LIMIT(int64)
DEFINE_INTTYPE_LIMIT(uint64)

#undef DEFINE_INTTYPE_LIMIT

} // namespace std

#endif
