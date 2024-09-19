// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/contrib/spinlessantifermion-u1-schwinger.h
//
// Copyright (C) 2022 Ian McCulloch <ian@qusim.net>
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
#include "quantumnumbers/u1.h"

inline
LatticeSite SpinlessAntifermionU1(std::string const& Sym1 = "N",
                              std::string const& ParityOp = "P")
{
   SymmetryList Symmetry(Sym1+":U(1)");
   QuantumNumbers::QNConstructor<QuantumNumbers::U1> QN(Symmetry);
   SiteBasis Basis(Symmetry);
   SiteOperator CH, C, P, R, N, Nf, Np, I;
   LatticeSite Site("U(1) Antifermion");

   // Following the Dirac Sea interpretation, the vacuum state has antifermion
   // sites fully occupied, so the antifermion is a 'hole', meaning the absence of a particle.
   // To keep the vacuum state as quantum number zero however we do a shift so that an occuped
   // sea state has quantum number 0.
   Basis.push_back("sea", QN(0));
   Basis.push_back("hole", QN(-1));

   OperatorDescriptions OpDescriptions;
   OpDescriptions.add_operators()
      ("I"       , "identity")
      ("R"       , "reflection")
      ("P"       , "fermion parity")
      ("C"       , "create antifermion")
      ("CH"      , "annihilate antifermion")
      ("N"       , "number operator = CH*C (1 for a sea state, 0 for a hole)")
      ("Nf"      , "number of fermions (-1 for an antifermion)")
      ("Np"      , "number of particles")
      ;

   LatticeCommute Fermionic = ParityOp == "P" ? LatticeCommute::Fermionic : LatticeCommute(ParityOp);

   CH = SiteOperator(Basis, QN(1), Fermionic);
   C = SiteOperator(Basis, QN(-1), Fermionic);
   P = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   R = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   N = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   Nf = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   Np = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   I = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);

   // annihilate antifermion, or alternatively: construct a sea particle from a hole
   CH("sea", "hole") = 1;

   // create antifermion
   C = adjoint(CH);

   // parity = (-1)^N - this is actually ambiguous, but the overall sign doesn't matter anyway
   P("sea"  , "sea")  =  1;
   P("hole" , "hole") = -1;

   // spatial reflection
   R("sea"  , "sea")  = 1;
   R("hole" , "hole") = 1;

   // particle number = CH*C.  This is +1 for a sea state (vacuum) and 0 for a hole (antifermion)
   N = CH*C;

   // Number of fermions.  This is -1 for a hole (antiparticle)
   Nf("hole", "hole") = -1;

   // Number of particles = +1 for both a fermion and an antifermion
   Np("hole", "hole") = 1;

   // identity
   I("sea",  "sea")  = 1;
   I("hole", "hole") = 1;

   Site["I"] = I;
   Site[ParityOp] = P;
   Site["N"] = N;
   Site["Nf"] = Nf;
   Site["Np"] = Np;
   Site["R"] = R;
   Site["CH"] = CH;
   Site["C"] = C;

   Site.set_operator_descriptions(OpDescriptions);

   return Site;
}
