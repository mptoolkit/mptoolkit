// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/clock-zn.h
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

#if !defined(MPTOOLKIT_MODELS_CONTRIB_CLOCK_ZN_H)
#define MPTOOLKIT_MODELS_CONTRIB_CLOCK_ZN_H

#include "common/stringutil.h"
#include "lattice/latticesite.h"
#include <complex>
#include <stdexcept>
#include <string>

inline
QuantumNumbers::QuantumNumber
ZnClockQuantumNumber(SymmetryList const& Symmetry, int q, int charge)
{
   int CanonicalCharge = ((charge % q) + q) % q;
   return QuantumNumbers::QuantumNumber(Symmetry, ConvertToString(CanonicalCharge));
}

inline
LatticeSite ZnClockSite(int q)
{
   if (q < 2)
      throw std::runtime_error("Z_q clock site supports q >= 2");

   SymmetryList Symmetry("Q:Z_" + ConvertToString(q));
   SiteBasis Basis(Symmetry);
   LatticeSite Site("Z_" + ConvertToString(q) + " clock");

   for (int s = 0; s < q; ++s)
   {
      Basis.push_back(std::to_string(s), ZnClockQuantumNumber(Symmetry, q, s));
   }

   OperatorDescriptions OpDescriptions;
   // Fourier-transformed relative to the plain Potts site: Gamma/X is diagonal
   // and neutral, while Omega/Z shifts between Z_q charge sectors.
   OpDescriptions.add_operators()
      ("I"     , "identity")
      ("R"     , "reflection")
      ("P"     , "fermion parity")
      ("Omega" , "charge +1 shift operator in the Gamma-diagonal basis (alias Z)")
      ("OmegaD", "charge -1 adjoint of Omega (alias ZD)")
      ("Gamma" , "neutral diagonal symmetry generator in this basis (alias X)")
      ("GammaD", "neutral adjoint of Gamma (alias XD)")
      ("Z"     , "alias for Omega")
      ("ZD"    , "alias for OmegaD")
      ("X"     , "alias for Gamma")
      ("XD"    , "alias for GammaD")
      ;

   QuantumNumbers::QuantumNumber q0 = ZnClockQuantumNumber(Symmetry, q, 0);
   QuantumNumbers::QuantumNumber q1 = ZnClockQuantumNumber(Symmetry, q, 1);
   QuantumNumbers::QuantumNumber qm1 = ZnClockQuantumNumber(Symmetry, q, -1);

   SiteOperator I(Basis, q0, LatticeCommute::Bosonic);
   SiteOperator P(Basis, q0, LatticeCommute::Bosonic);
   SiteOperator R(Basis, q0, LatticeCommute::Bosonic);
   SiteOperator Omega(Basis, q1, LatticeCommute::Bosonic);
   SiteOperator OmegaD(Basis, qm1, LatticeCommute::Bosonic);
   SiteOperator Gamma(Basis, q0, LatticeCommute::Bosonic);
   SiteOperator GammaD(Basis, q0, LatticeCommute::Bosonic);

   for (int s = 0; s < q; ++s)
   {
      std::string ss = std::to_string(s);
      std::string sp = std::to_string((s+1)%q);
      std::string sm = std::to_string((s+q-1)%q);
      std::complex<double> phase = std::exp(std::complex<double>(0.0, 1.0) *
                                            2.0 * math_const::pi * double(s) / double(q));

      I(ss, ss) = 1.0;
      P(ss, ss) = 1.0;
      R(ss, ss) = 1.0;
      Omega(sp, ss) = 1.0;
      OmegaD(sm, ss) = 1.0;
      Gamma(ss, ss) = phase;
      GammaD(ss, ss) = std::conj(phase);
   }

   Site.arg("q") = q;

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

#endif
