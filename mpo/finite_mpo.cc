// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mpo/finite_mpo.cc
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

inline
std::ostream&
operator<<(std::ostream& out, FiniteMPO const& x)
{
   return out << x.data();
}

inline
FiniteMPO::FiniteMPO(GenericMPO const& Other)
   : Data(Other)
{
   //   CHECK(Data.back().Basis2().is_identity())("Finite operator: right basis must be scalar");
   CHECK(Data.back().Basis2().is_regular())("Finite operator: left basis must be regular");
}
