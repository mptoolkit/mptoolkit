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
   SiteOperator Omega, OmegaD, Gamma, GammaD, R, P, I;
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
      ("Omega", "clock operator, diagonal in the Potts basis (alias Z)")
      ("OmegaD", "adjoint of Omega (alias ZD)")
      ("Gamma", "shift operator (alias X)")
      ("GammaD", "adjoint of Gamma (alias XD)")
      ("Z"   , "alias for Omega")
      ("ZD"  , "alias for OmegaD")
      ("X"   , "alias for Gamma")
      ("XD"  , "alias for GammaD")
      ;

   P = SiteOperator(Basis, q, LatticeCommute::Bosonic);
   R = SiteOperator(Basis, q, LatticeCommute::Bosonic);
   I = SiteOperator(Basis, q, LatticeCommute::Bosonic);
   Omega = SiteOperator(Basis, q, LatticeCommute::Bosonic);
   OmegaD = SiteOperator(Basis, q, LatticeCommute::Bosonic);
   Gamma = SiteOperator(Basis, q, LatticeCommute::Bosonic);
   GammaD = SiteOperator(Basis, q, LatticeCommute::Bosonic);

   for (int s = 0; s < qq; ++s)
   {
      std::string ss = std::to_string(s);
      std::string sp = std::to_string((s+1)%qq);
      std::string sm = std::to_string((s+qq-1)%qq);
      std::complex<double> phase = std::exp(std::complex<double>(0.0, 1.0) *
                                            2.0 * math_const::pi * double(s) / double(qq));
      I(ss, ss) = 1.0;
      P(ss, ss) = 1.0;
      R(ss, ss) = 1.0;
      Omega(ss, ss) = phase;
      OmegaD(ss, ss) = std::conj(phase);
      Gamma(sp, ss) = 1;
      GammaD(sm, ss) = 1;
   }

   Site.arg("q") = qq;

   Site["I"] = I;
   Site["P"] = P;
   Site["R"] = R;
   Site["Omega"] = Omega;
   Site["OmegaD"] = OmegaD;
   Site["Gamma"] = Gamma;
   Site["GammaD"] = GammaD;
   Site["Z"] = Omega;
   Site["ZD"] = OmegaD;
   Site["X"] = Gamma;
   Site["XD"] = GammaD;

   Site.set_operator_descriptions(OpDescriptions);

   return Site;
}
