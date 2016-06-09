// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/vector_utility.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#if !defined(VECTOR_UTILITY_H_JSHDCUIH78943QY7YPO7YP89)
#define VECTOR_UTILITY_H_JSHDCUIH78943QY7YPO7YP89

#include "vector.h"
#include "matrix_utility.h" // to get the random functions
#include <stdlib.h>
#include <complex>

namespace LinearAlgebra
{

template <typename Func>
Vector<typename Func::value_type>
generate_vector(size_type Size, Func f = Func())
{
   Vector<typename Func::value_type> Result(Size);
   for (size_type i = 0; i < Size; ++i)
   {
      Result[i] = f();
   }
   return Result;
}

template <typename Scalar>
Vector<Scalar>
generate_vector(size_type Size, Scalar (&f)())
{
   Vector<typename make_value<Scalar>::type> Result(Size);
   for (size_type i = 0; i < Size; ++i)
   {
      Result[i] = f();
   }
   return Result;
}

template <typename Scalar>
Vector<Scalar> random_vector(size_type Size)
{
   return generate_vector(Size, random<Scalar>);
}

} // namespace LinearAlgebra

#endif
