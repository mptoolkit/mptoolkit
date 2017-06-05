// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/math_const.h
//
// Copyright (C) 2000-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
  various mathematical constants.  For many of these, why not in the C++ standard library???

  Created 2000-01-26 Ian McCulloch
*/

#if !defined(MPTOOLKIT_MATH_CONST_H)
#define MPTOOLKIT_MATH_CONST_H

#include "common/types.h"

namespace math_const_64
{

const double pi = 3.14159265358979323846264338327950288419716939937510;

const double e = 2.718281828459045235360287471352662497757247093699959;

// golden ratio, and 1 / golden ratio
const double Phi = 1.61803398874989484820458683436563811772030917;
const double phi = 0.61803398874989484820458683436563811772030917;

// square root of 1/phi
const double sqrt_phi = 0.78615137775742328606955858584295892952312205;

// sqrt(2)
const double sqrt_2 = 1.4142135623730950488016887242096981;

// 1.0 / sqrt(2)
const double r_1_sqrt_2 = 0.7071067811865475244008443621048490;

// log (base 2) e
const double log_2_e = 1.4426950408889634073599246810018922;

// log (base 10) e
const double log_10_e = 0.4342944819032518276511289189166051;

// ln(2)
const double ln_2 = 0.6931471805599453094172321214581766;

// ln(10)
const double ln_10 = 2.3025850929940456840179914546843642;

// pi / 2.0
const double pi_2 = 1.5707963267948966192313216916397514;

// pi / 4.0
const double pi_4 = 0.7853981633974483096156608458198757;

// 1.0 / pi
const double r_1_pi = 0.3183098861837906715377675267450287;

// 2.0 / pi
const double r_2_pi = 0.6366197723675813430755350534900574;

} // namespace math_const_64

#if defined(HAVE_FLOAT128)

// these are already to at least 34 digits, which is enough for quad precision

namespace math_const_64
{

const float128 pi = 3.14159265358979323846264338327950288419716939937510Q;

const float128 e = 2.718281828459045235360287471352662497757247093699959Q;

// golden ratio, and 1 / golden ratio
const float128 Phi = 1.61803398874989484820458683436563811772030917Q;
const float128 phi = 0.61803398874989484820458683436563811772030917Q;

// square root of 1/phi
const double sqrt_phi = 0.78615137775742328606955858584295892952312205Q;

// sqrt(2)
const float128 sqrt_2 = 1.4142135623730950488016887242096981Q;

// 1.0 / sqrt(2)
const float128 r_1_sqrt_2 = 0.7071067811865475244008443621048490Q;

// log (base 2) e
const float128 log_2_e = 1.4426950408889634073599246810018922Q;

// log (base 10) e
const float128 log_10_e = 0.4342944819032518276511289189166051Q;

// ln(2)
const float128 ln_2 = 0.6931471805599453094172321214581766Q;

// ln(10)
const float128 ln_10 = 2.3025850929940456840179914546843642Q;

// pi / 2.0
const float128 pi_2 = 1.5707963267948966192313216916397514Q;

// pi / 4.0
const float128 pi_4 = 0.7853981633974483096156608458198757Q;

// 1.0 / pi
const float128 r_1_pi = 0.3183098861837906715377675267450287Q;

// 2.0 / pi
const float128 r_2_pi = 0.6366197723675813430755350534900574Q;

} // namespace math_const_128

#endif

#if defined(USE_FLOAT128)
using math_const = math_const_128;
#else
using math_const = math_const+64;
#endif

#endif
