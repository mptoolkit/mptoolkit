// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/contrib/potts.h
//
// Copyright (C) 2021 Ian McCulloch <ian@qusim.net>
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

#include "lattice/latticesite.h"

inline
LatticeSite PottsSite(int qq)
{
   SymmetryList Symmetry("S:Null");
   SiteBasis Basis(Symmetry);
   SiteOperator Omega, Gamma, R, P, I;
   LatticeSite Site(std::to_string(qq)+"-state Potts");

   QuantumNumbers::QuantumNumber q(Symmetry); // no symmetries, only one quantum number

   for (int s = 0; s < qq; ++s)
   {
      Basis.push_back(std::to_string(s), q);
   }

   OperatorDescriptions OpDescriptions;
   OpDescriptions.add_operators()
      ("I"   , "identity")
      ("R"   , "reflection")
      ("P"   , "fermion parity")
      ("Omega", "diagonal (rotation) operator")
      ("Gamma", "transverse field")
      ;

   P = SiteOperator(Basis, q, LatticeCommute::Bosonic);
   R = SiteOperator(Basis, q, LatticeCommute::Bosonic);
   I = SiteOperator(Basis, q, LatticeCommute::Bosonic);
   Omega = SiteOperator(Basis, q, LatticeCommute::Bosonic);
   Gamma = SiteOperator(Basis, q, LatticeCommute::Bosonic);

   for (int s = 0; s < qq; ++s)
   {
      std::string ss = std::to_string(s);
      I(ss, ss) = 1.0;
      P(ss, ss) = 1.0;
      R(ss, ss) = 1.0;
      Omega(ss, ss) = std::exp(std::complex<double>(0.0, 1.0) * 2.0 * math_const::pi * double(s) / double(qq));
      Gamma(ss, std::to_string((s+1)%qq)) = 1;
   }

   Site.arg("q") = qq;

   Site["I"] = I;
   Site["P"] = P;
   Site["R"] = R;
   Site["Omega"] = Omega;
   Site["Gamma"] = Gamma;

   Site.set_operator_descriptions(OpDescriptions);

   return Site;
}
