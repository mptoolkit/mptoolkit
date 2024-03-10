// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// common/math_const.h
//
// Copyright (C) 2000-2016 Ian McCulloch <ian@qusim.net>
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

/*
  various mathematical constants.  For many of these, why not in the C++ standard library???

  Created 2000-01-26 Ian McCulloch
*/

#if !defined(MATH_CONST_H_T7845TY7ERY7439JYQWTYJJKHF34Y789W3YUSDHE78)
#define MATH_CONST_H_T7845TY7ERY7439JYQWTYJJKHF34Y789W3YUSDHE78

namespace math_const
{

// Just in case this is ever compiled on a machine with HUGE floating point,
// we've got 100 decimal places :)
const double pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679;

// We've got e to 100 places too
const double e = 2.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274;

// golden ratio, and 1 / golden ratio
const double Phi = 1.61803398874989484820458683436563811772030917980576286213544862270526046281890;
const double phi = 0.61803398874989484820458683436563811772030917980576286213544862270526046281890;

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

} // namespace math_const

#endif
