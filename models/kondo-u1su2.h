// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/kondo-u1su2.h
//
// Copyright (C) 2012-2026 Ian McCulloch <ian@qusim.net>
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
// Site block for a hubbard-like site with an impurity spin

#ifndef MPTOOLKIT_MODELS_KONDO_U1SU2_H
#define MPTOOLKIT_MODELS_KONDO_U1SU2_H

#include "models/kondo-common.h"
#include "models/fermion-u1su2.h"

inline LatticeSite KondoU1SU2(half_int Spin, std::string const& Sym1 = "N",
                            std::string const& Sym2 = "S")
{
   LatticeSite Electron = FermionU1SU2(Sym1, Sym2);
   // The Kondo site has always called the number operator N, even with a
   // custom charge-symmetry name.
   return KondoDetail::MakeSite(Electron, Spin, Sym2, "U(1)xSU(2) Kondo impurity", Sym1);
}

// Preserve source compatibility as well as the default spin-1/2 site.
inline LatticeSite KondoU1SU2(std::string const& Sym1 = "N", std::string const& Sym2 = "S")
{
   return KondoU1SU2(half_int(0.5), Sym1, Sym2);
}

#endif
