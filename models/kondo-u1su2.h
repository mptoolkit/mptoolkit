// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/kondo-u1su2.h
//
// Copyright (C) 2012-2017 Ian McCulloch <ian@qusim.net>
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

#include "lattice/latticesite.h"
#include "quantumnumbers/u1.h"
#include "quantumnumbers/su2.h"

inline
LatticeSite
KondoU1SU2(std::string const& Sym1 = "N", std::string const& Sym2 = "S")
{
   SymmetryList Symmetry(Sym1+":U(1),"+Sym2+":SU(2)");
   QuantumNumbers::QNConstructor<QuantumNumbers::U1,QuantumNumbers::SU2> QN(Symmetry);
   SiteBasis Basis(Symmetry);
   SiteOperator C, CH, P, R, N, S, Sf, Sc, ScSf, I, Hu, Pdouble,
      Ep, Em, Ns, Nh, Pg, CP, CHP;
   LatticeSite Site("U(1)xSU(2) Kondo impurity");

   Basis.push_back("empty",  QN(0, 0.5));
   Basis.push_back("double", QN(2, 0.5));
   Basis.push_back("singlet", QN(1, 0));
   Basis.push_back("triplet", QN(1, 1));

   OperatorDescriptions OpDescriptions;
   OpDescriptions.add_operators()
      ("I"       , "identity")
      ("R"       , "reflection")
      ("C"       , "annihilation operator")
      ("CH"      , "creation operator")
      ("P"       , "fermion parity")
      ("N"       , "number operator")
      ("N_S"     , "number of spins")
      ("N_H"     , "number of holons")
      ("S"       , "spin vector")
      ("Sc"      , "conduction spin vector")
      ("Sf"      , "f-spin vector") 
      ("Ep"      , "create double-occupied site (eta paring)")
      ("Em"      , "annihilate double-occupied site (eta paring)")
      ("Hu"      , "symmetrized Coulomb operator (n_up - 1/2) * (n_down - 1/2)")
      ("Pdouble" , "projector onto the double-occupied site")
      ("ScSf"    , "Kondo spin interaction Sc . Sf")
      ("Pg"      , "Gutswiller projector = 1-Pdouble")
      ;

   C = SiteOperator(Basis, QN(-1,0.5), LatticeCommute::Fermionic);
   CH = SiteOperator(Basis, QN(1,0.5), LatticeCommute::Fermionic);
   P = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   R = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   N = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   Hu = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   Pdouble = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   Pg = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   S = SiteOperator(Basis, QN(0,1), LatticeCommute::Bosonic);
   Sf = SiteOperator(Basis, QN(0,1), LatticeCommute::Bosonic);
   Sc = SiteOperator(Basis, QN(0,1), LatticeCommute::Bosonic);
   ScSf = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   I = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   Ns = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   Nh = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);

   // annihilate fermion
   C("empty",  "singlet")    =  sqrt(0.5);
   C("empty",  "triplet")    =  sqrt(1.5);
   C("singlet", "double")    = -1.0;
   C("triplet", "double")    =  1.0;

   // create fermion
   CH = adjoint(C);

   // create double occupied site (eta pairing)
   Ep = sqrt(2.0) * prod(CH, CH, QN(2,0));
   Em = sqrt(2.0) * prod(C, C, QN(-2,0));

   // parity = (-1)^N
   P("empty",   "empty")      =  1;
   P("singlet", "singlet")    = -1;
   P("triplet", "triplet")    = -1;
   P("double",  "double")     =  1;

   // spatial reflection
   R("empty",  "empty")       =  1;
   R("singlet", "singlet")    =  1;
   R("triplet", "triplet")    =  1;
   R("double", "double")      = -1;

   // particle number
   N("singlet", "singlet")    =  1;
   N("triplet", "triplet")    =  1;
   N("double",  "double")     =  2;

   // symmetrized coulomb operator = (n_up - 1/2) * (n_down - 1/2)
   Hu("empty",  "empty")       =  0.25;
   Hu("singlet", "singlet")    = -0.25;
   Hu("triplet", "triplet")    = -0.25;
   Hu("double", "double")      =  0.25;

   // projection onto double occupied state == asymmetric coulomb operator == n_up * n_down
   Pdouble("double", "double") = 1.0;

   // identity
   I("empty",  "empty")       =  1;
   I("singlet", "singlet")    =  1;
   I("triplet", "triplet")    =  1;
   I("double", "double")      =  1;

   Pg = I - Pdouble; // Gutzwiller projector

   // S
   S("empty",   "empty")     = std::sqrt(0.75);
   S("double",  "double")    = std::sqrt(0.75);
   S("triplet", "triplet")   = std::sqrt(2.0);

   // conduction spin
   Sc("singlet", "triplet")   = -std::sqrt(0.75);
   Sc("triplet", "singlet")   =  0.5;
   Sc("triplet", "triplet")   =  std::sqrt(0.5);

   // local spin
   Sf("empty",   "empty")     =  std::sqrt(0.75);
   Sf("double",  "double")    =  std::sqrt(0.75);

   Sf("singlet", "triplet")   =  std::sqrt(0.75);
   Sf("triplet", "singlet")   = -0.5;
   Sf("triplet", "triplet")   = std::sqrt(0.5);

   // number of spins in the conduction band
   Ns("singlet", "singlet") = 1;
   Ns("triplet", "triplet") = 1;

   // Local spin-spin interaction, equivalent to -sqrt(3.0) * prod(Sc, Sf)
   ScSf("singlet", "singlet") = -0.75;
   ScSf("triplet", "triplet") = 0.25;
   // number of holons
   Nh = I - Ns;

   Site["I"] = I;
   Site["P"] = P;
   Site["N"] = N;
   Site["Hu"] = Hu;
   Site["Pdouble"] = Pdouble;
   Site["Pg"] = Pg;
   Site["S"] = S;
   Site["Sf"] = Sf;
   Site["Sc"] = Sc;
   Site["C"] = C;
   Site["R"] = R;
   Site["CH"] = CH;
   Site["Ep"] = Ep;
   Site["Em"] = Em;
   Site["N_S"] = Ns;
   Site["N_H"] = Nh;
   Site["ScSf"] = ScSf;

   Site.set_operator_descriptions(OpDescriptions);

   return Site;
}
