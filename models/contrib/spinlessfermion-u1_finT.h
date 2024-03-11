// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/contrib/spinlessfermion-u1_finT.h
//
// Copyright (C) 2012-2021 Ian McCulloch <ian@qusim.net>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Reseach publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

#include "lattice/latticesite.h"
#include "quantumnumbers/u1.h"

inline
LatticeSite SpinlessFermionU1(std::string const& Sym1 = "N",
                              std::string const& ParityOp = "P")
{
   SymmetryList Symmetry(Sym1+":U(1)");
   QuantumNumbers::QNConstructor<QuantumNumbers::U1> QN(Symmetry);
   SiteBasis Basis(Symmetry);
   SiteOperator CH, C, P, R, N, I, P_0, P_1;
   LatticeSite Site("U(1) Fermion");

   Basis.push_back("empty", QN(0));
   Basis.push_back("single", QN(1));

   OperatorDescriptions OpDescriptions;
   OpDescriptions.add_operators()
      ("I"       , "identity")
      ("R"       , "reflection")
      ("P"       , "fermion parity")
      ("C"       , "annihilation fermion")
      ("CH"      , "creation fermion")
      ("N"       , "number operator")
      ("P_0"     , "Project onto empty site")
      ("P_1"     , "Project onto singly-occupied site")
      ;

   LatticeCommute Fermionic = ParityOp == "P" ? LatticeCommute::Fermionic : LatticeCommute(ParityOp);

   CH = SiteOperator(Basis, QN(1), Fermionic);
   C = SiteOperator(Basis, QN(-1), Fermionic);
   P = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   R = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   N = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   I = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   P_0 = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   P_1 = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);

   // annihilate fermion
   C("empty", "single") = 1;

   // create fermion
   CH = adjoint(C);

   // parity = (-1)^N
   P("empty",  "empty")  =  1;
   P("single", "single") = -1;

   // spatial reflection
   R("empty",  "empty")  = 1;
   R("single", "single") = 1;

   // particle number
   N("single", "single") = 1;

   // identity
   I("empty",  "empty")  = 1;
   I("single", "single") = 1;

   P_0("empty", "empty") = 1;
   P_1("single", "single") = 1;

   Site["I"] = I;
   Site[ParityOp] = P;
   Site["N"] = N;
   Site["R"] = R;
   Site["CH"] = CH;
   Site["C"] = C;
   Site["P_0"] = P_0;
   Site["P_1"] = P_1;

   Site.set_operator_descriptions(OpDescriptions);

   return Site;
}
