// -*- C++ -*- $Id$

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
