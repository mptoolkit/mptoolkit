// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/kondo-su2.h
//
// Copyright (C) 2026 Ian McCulloch <ian@qusim.net>
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

#ifndef MPTOOLKIT_MODELS_KONDO_SU2_H
#define MPTOOLKIT_MODELS_KONDO_SU2_H

#include "models/kondo-common.h"
#include "models/fermion-su2.h"

inline LatticeSite KondoSU2(half_int Spin = half_int(0.5), std::string const& Sym = "S")
{
   // Fermion parity is retained as an operator (and in the commutation
   // metadata), even though particle number is not a symmetry label.
   return KondoDetail::MakeSite(FermionSU2(Sym), Spin, Sym, "SU(2) Kondo impurity");
}

#endif
