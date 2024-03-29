// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mpo/basic_finite_mpo.cc
//
// Copyright (C) 2013-2017 Ian McCulloch <ian@qusim.net>
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

inline
std::ostream&
operator<<(std::ostream& out, BasicFiniteMPO const& x)
{
   return out << x.data();
}

inline
BasicFiniteMPO::BasicFiniteMPO(GenericMPO const& Other)
   : Data(Other)
{
   //   CHECK(Data.back().Basis2().is_identity())("Finite operator: right basis must be scalar");
   CHECK(Data.back().Basis2().is_regular())("Finite operator: left basis must be regular");
}

inline
BasicFiniteMPO coarse_grain(BasicFiniteMPO const& Op, int N)
{
   return BasicFiniteMPO(coarse_grain(Op.data(), N));
}
