// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// blas/functors.h
//
// Copyright (C) 2017 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

// Generic function objects for wrapping operators and functions.
// These are useful for nested operations, eg a nested tensor product can
// be evaluated as tensor_prod(A, B, TensorProd())
// versus the default tensor_prod(A, B, Multiplication())

#if !defined(MPTOOLKIT_BLAS_FUNCTORS_H)
#define MPTOOLKIT_BLAS_FUNCTORS_H

namespace blas
{

struct Multiplication
{
   template <typename T, typename U>
   auto operator()(T const& x, U const& y) const
   {
      return x*y;
   }
};

struct Addition
{
   template <typename T, typename U>
   auto operator()(T const& x, U const& y) const
   {
      return x+y;
   }
};

struct Subtraction
{
   template <typename T, typename U>
   auto operator()(T const& x, U const& y) const
   {
      return x-y;
   }
};

} // namespace blas

#endif
