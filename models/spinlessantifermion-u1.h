// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/spinlessantifermion-u1.h
//
// Copyright (C) 2012-2016 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2022 Jesse Osborne <j.osborne@uqconnect.edu.au>
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
                                  std::string const& ParityOp = "P",
                                  bool Bosonic = false)
{
   SymmetryList Symmetry(Sym1+":U(1)");
   QuantumNumbers::QNConstructor<QuantumNumbers::U1> QN(Symmetry);
   SiteBasis Basis(Symmetry);
   SiteOperator CH, C, P, R, N, I;
   LatticeSite Site("U(1) Antifermion");

   Basis.push_back("empty", QN(0));
   Basis.push_back("single", QN(-1));

   OperatorDescriptions OpDescriptions;
   OpDescriptions.add_operators()
      ("I"       , "identity")
      ("R"       , "reflection")
      ("P"       , "antifermion parity")
      ("C"       , "creation antifermion")
      ("CH"      , "annihilation antifermion")
      ("N"       , "number operator")
      ;

   LatticeCommute Fermionic = ParityOp == "P" ? (Bosonic ? LatticeCommute::Bosonic : LatticeCommute::Fermionic)
                                              : LatticeCommute(ParityOp);

   CH = SiteOperator(Basis, QN(1), Fermionic);
   C = SiteOperator(Basis, QN(-1), Fermionic);
   P = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   R = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   N = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   I = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);

   // create antifermion
   C("single", "empty") = 1;

   // annihilate antifermion
   CH = adjoint(C);

   // parity = (-1)^N
   P("empty",  "empty")  =  1;
   P("single", "single") = -1;

   // spatial reflection
   R("empty",  "empty")  = 1;
   R("single", "single") = 1;

   // particle number
   N("single", "single") = -1;

   // identity
   I("empty",  "empty")  = 1;
   I("single", "single") = 1;

   Site["I"] = I;
   Site[ParityOp] = P;
   Site["N"] = N;
   Site["R"] = R;
   Site["CH"] = CH;
   Site["C"] = C;

   Site.set_operator_descriptions(OpDescriptions);

   return Site;
}
