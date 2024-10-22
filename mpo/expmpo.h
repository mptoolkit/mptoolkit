// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mpo/basic_triangular_mpo.h
//
// Copyright (C) 2023 Jesse Osborne <j.osborne@uqconnect.edu.au>
// Copyright (C) 2024 Ian McCulloch <ian@qusim.net>
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
// Exponential of an MPO, based on https://arxiv.org/abs/2302.14181

#if !defined(MPTOOLKIT_MPO_EXPMPO_H)
#define MPTOOLKIT_MPO_EXPMPO_H

#include "common/namedenum.h"
#include "basic_triangular_mpo.h"
#include "product_mpo.h"
#include <array>

// Optimized first-order scheme, calculates Result' = exp(Tau * x) to first order, including
// additional terms where possible, that don't increase the dimension of the MPO
ProductMPO expmpo_first(BasicTriangularMPO const& x, std::complex<double> Tau = 1.0, int DysonOrder = 6);

// Strict first-order scheme, calculates Result' = exp(Tau * x) truncating all expansions to first order
ProductMPO expmpo_first_strict(BasicTriangularMPO const& x, std::complex<double> Tau = 1.0);

// Optimized second-order scheme
ProductMPO expmpo_second(BasicTriangularMPO const& x, std::complex<double> Tau = 1.0, int DysonOrder = 6);

//
// Generic version using a named enumeration to choose the scheme
//

struct ExpMpoSchemeTraits
{
   enum Enum { first, first_strict, second };
   static constexpr Enum Default = first;
   static constexpr char const* StaticName = "expmpo scheme";
   static constexpr std::array<char const*, 3> Names = { "1", "1strict", "2" };
};

struct ExpMpoSchemeParameters
{
   ExpMpoSchemeParameters() : DysonOrder(6) {}
   int DysonOrder;
};

using ExpMpoScheme = NamedEnumeration<ExpMpoSchemeTraits>;

ProductMPO expmpo(BasicTriangularMPO const& x, ExpMpoScheme Scheme = ExpMpoScheme(), ExpMpoSchemeParameters p = ExpMpoSchemeParameters());

#endif
