// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/spin-su2.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "lattice/latticesite.h"
#include "quantumnumbers/su2.h"

inline
LatticeSite
SpinSU2(half_int Spin, std::string const& Sym = "S")
{
   SymmetryList Symmetry(Sym+":SU(2)");
   QuantumNumbers::QNConstructor<QuantumNumbers::SU2> QN(Symmetry);
   SiteBasis Basis(Symmetry);
   SiteOperator S, Sn, R, P, I, Q, T, F;
   LatticeSite Site("SU(2) Spin "+to_string_fraction(Spin));

   std::string SpinName = boost::lexical_cast<std::string>(Spin);
   Basis.push_back(SpinName, QN(Spin));

   OperatorDescriptions OpDescriptions;
   OpDescriptions.add_operators()
      ("I"   , "identity")
      ("R"   , "reflection")
      ("P"   , "fermion parity")
      ("S"   , "spin vector")
      ("Q"   , "spin tensor [spin 2]")
      ("T"   , "spin tensor [spin 3]")
      ("F"   , "spin tensor [spin 4]")
      ;


   S = SiteOperator(Basis, QN(1), LatticeCommute::Bosonic);
   P = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   R = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   I = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);

   S(SpinName, SpinName) = sqrt(Spin * (Spin+1));
   I(SpinName, SpinName) = 1.0;
   P = I;
   R = I;

   Q = prod(S, S, QN(2));
   T = prod(Q, S, QN(3));
   F = prod(Q, Q, QN(4));
   
   Site.arg("Spin") = Spin.to_double();

   Site["I"] = I;
   Site["P"] = P;
   Site["R"] = R;
   Site["S"] = S;
   Site["Q"] = Q;
   Site["T"] = T;
   Site["F"] = F;

   Site.set_operator_descriptions(OpDescriptions);

   return Site;
}
