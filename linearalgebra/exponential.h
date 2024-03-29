// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/exponential.h
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
//
// Exponential function for matrices.

#if !defined(EXPONENTIAL_H_SDJKCH48935Y78Y7RFEH89P)
#define EXPONENTIAL_H_SDJKCH48935Y78Y7RFEH89P

#include "matrix.h"

namespace LinearAlgebra
{


//
// Exponentiate
//
// Calculate exp(M) of a matrix M, using Expokit.
//

template <typename M, typename Mi = typename interface<M>::type>
struct ImplementExponentiate {};

template <typename M>
inline
typename ImplementExponentiate<M>::result_type
Exponentiate(double t, M const& m)
{
   return ImplementExponentiate<M>()(t, m);
}

} // namespace LinearAlgebra

#include "exponential.cc"

#endif
