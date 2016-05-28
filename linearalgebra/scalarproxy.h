// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/scalarproxy.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

/*
  Created 2004-06-22 Ian McCulloch

  Defines function scalar(T) to create a proxy to treat T as a scalar.

  Usage:
  Matrix<double> M;
  std::complex<double> x;
  scalar(x) * M;   // equivalent to left_scalar_prod(x, M)
*/

#if !defined(SCALARPROXY_H_DSHFIUYR984Y987HYIUAH874WLO)
#define SCALARPROXY_H_DSHFIUYR984Y987HYIUAH874WLO

namespace LinearAlgebra
{

// The actual proxy class goes in namespace Private -
// the referent object is held by reference (so goes out of scope at the end of the
// expression where it is defined), so we don't want to use these objects outside of
// expressions.
namespace Private
{

template <typename T>
struct ScalarProxy
{
   typedef T const value_type;
   ScalarProxy(T const& value_) : value(value_) {}

   operator T const&() const { return value; }

   T const& value;
};

} // namespace Private

template <typename T>
inline
Private::ScalarProxy<T>
scalar(T const& x)
{
   return Private::ScalarProxy<T>(x);
}

} // namespace LinearAlgebra

#endif
