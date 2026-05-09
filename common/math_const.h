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

#include <numbers>

namespace math_const
{

inline constexpr double pi = std::numbers::pi_v<double>;

inline constexpr double e = std::numbers::e_v<double>;

// golden ratio, and 1 / golden ratio
inline constexpr double Phi = std::numbers::phi_v<double>;
inline constexpr double phi = 1.0 / Phi;

// sqrt(2)
inline constexpr double sqrt_2 = std::numbers::sqrt2_v<double>;

// 1.0 / sqrt(2)
inline constexpr double r_1_sqrt_2 = 1.0 / sqrt_2;

// log (base 2) e
inline constexpr double log_2_e = std::numbers::log2e_v<double>;

// log (base 10) e
inline constexpr double log_10_e = std::numbers::log10e_v<double>;

// ln(2)
inline constexpr double ln_2 = std::numbers::ln2_v<double>;

// ln(10)
inline constexpr double ln_10 = std::numbers::ln10_v<double>;

// pi / 2.0
inline constexpr double pi_2 = pi / 2.0;

// pi / 4.0
inline constexpr double pi_4 = pi / 4.0;

// 1.0 / pi
inline constexpr double r_1_pi = std::numbers::inv_pi_v<double>;

// 2.0 / pi
inline constexpr double r_2_pi = 2.0 * r_1_pi;

} // namespace math_const

#endif
