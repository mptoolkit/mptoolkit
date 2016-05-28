// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/gmpalgebra.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER
//
// type traits for linear algebra library with GMP

#if !defined(GMPALGEBRA_H_HDS3247Y63478YFUIO)
#define GMPALGEBRA_H_HDS3247Y63478YFUIO

#include "common/gmpint.h"
#include "common/gmprational.h"
#include "scalar.h"

namespace LinearAlgebra
{

template <>
struct PromoteTraits<gmp::bigint, gmp::bigint>
{
   typedef gmp::bigint first_argument_type;
   typedef gmp::bigint second_argument_type;
   typedef gmp::bigint value_type;
   typedef gmp::bigint result_type;
};

template <>
struct PromoteTraits<gmp::rational, gmp::rational>
{
   typedef gmp::rational first_argument_type;
   typedef gmp::rational second_argument_type;
   typedef gmp::rational value_type;
   typedef gmp::rational result_type;
};

// 1-norm and squared 2-norm

template <>
struct Norm2Sq<gmp::bigint, AnyScalar<gmp::bigint> >
{
   typedef gmp::bigint result_type;
   typedef gmp::bigint argument_type;
   result_type operator()(gmp::bigint const& x) const
   {
      return x*x;
   }
};

template <>
struct Norm2Sq<gmp::rational, AnyScalar<gmp::rational> >
{
   typedef gmp::rational result_type;
   typedef gmp::rational argument_type;
   result_type operator()(gmp::rational const& x) const
   {
      return x*x;
   }
};

template <>
struct NormFrobSq<gmp::bigint, AnyScalar<gmp::bigint> >
{
   typedef gmp::bigint result_type;
   typedef gmp::bigint argument_type;
   result_type operator()(gmp::bigint const& x) const
   {
      return x*x;
   }
};

template <>
struct NormFrobSq<gmp::rational, AnyScalar<gmp::rational> >
{
   typedef gmp::rational result_type;
   typedef gmp::rational argument_type;
   result_type operator()(gmp::rational const& x) const
   {
      return x*x;
   }
};

template <>
struct Norm1<gmp::bigint, AnyScalar<gmp::bigint> >
{
   typedef gmp::bigint result_type;
   typedef gmp::bigint argument_type;
   result_type operator()(gmp::bigint const& x) const
   {
      return abs(x);
   }
};

template <>
struct Norm1<gmp::rational, AnyScalar<gmp::rational> >
{
   typedef gmp::rational result_type;
   typedef gmp::rational argument_type;
   result_type operator()(gmp::rational const& x) const
   {
      return abs(x);
   }
};

// conjugation

template<>
struct ConjInterface<gmp::bigint, AnyScalar<gmp::bigint> >
   : Identity<gmp::bigint> {};

template<>
struct ConjInterface<gmp::rational, AnyScalar<gmp::rational> >
   : Identity<gmp::rational> {};

} // namespace LinearAlgebra

#endif
