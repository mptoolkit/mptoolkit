// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/old/kondo-so4.h
//
// Copyright (C) 2012-2016 Ian McCulloch <ian@qusim.net>
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
//
// Site block for a hubbard-like site with an impurity spin

#include "lattice/latticesite.h"
#include "quantumnumbers/u1.h"
#include "quantumnumbers/su2.h"




inline
LatticeSite CreateSO4KondoSiteCommon(std::string const& Sym1 = "Q", std::string const& Sym2 = "S")
{
   SymmetryList Symmetry(Sym1+":SU(2),"+Sym2+":SU(2)");
   QuantumNumbers::QNConstructor<QuantumNumbers::SU2,QuantumNumbers::SU2> QN(Symmetry);
   SiteBasis Basis(Symmetry);
   SiteOperator C_A, C_B, CH_A, CH_B, P, R, N, S, Sf, Sc, ScSf, I, Hu, Pdouble,
      Ep, Em, Ns, Nh, Pg;
   LatticeSite Site;

   Basis.push_back("singlet", QN(0, 0));
   Basis.push_back("triplet", QN(0, 1));
   Basis.push_back("holon", QN(0.5, 0.5));

   C_A = SiteOperator(Basis, QN(0.5,0.5), LatticeCommute::Fermionic);
   C_B = SiteOperator(Basis, QN(0.5,0.5), LatticeCommute::Fermionic);
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
   C_A("singlet", "holon")  =  sqrt(2.0);
   C_A("triplet", "holon")  =  sqrt(2.0);
   C_A("holon",  "singlet") = -sqrt(0.5);
   C_A("holon",  "triplet") =  sqrt(1.5);

   // annihilate fermion
   C_B("singlet", "holon")  =  sqrt(2.0);
   C_B("triplet", "holon")  =  sqrt(2.0);
   C_B("holon",  "singlet") =  sqrt(0.5);
   C_B("holon",  "triplet") = -sqrt(1.5);

   // create fermion
   CH_A = adjoint(C_A);
   CH_B = adjoint(C_B);

   // parity = (-1)^N
   P("singlet", "singlet")    = -1;
   P("triplet", "triplet")    = -1;
   P("holon",   "holon")      =  1;

   // parity modified creation/annihilation
   SiteOperator CP_A = C_A * P;
   SiteOperator CP_B = C_B * P;
   SiteOperator CHP_A = CH_A * P;
   SiteOperator CHP_B = CH_B * P;

   // symmetrized coulomb operator = (n_up - 1/2) * (n_down - 1/2)
   Hu("singlet", "singlet")    = -0.25;
   Hu("triplet", "triplet")    = -0.25;
   Hu("holon",   "holon")      =  0.25;

   // identity
   I("singlet", "singlet")    =  1;
   I("triplet", "triplet")    =  1;
   I("holon",  "holon")       =  1;

   // S
   S("triplet", "triplet")   = std::sqrt(2.0);
   S("holon",   "holon")     = std::sqrt(0.75);

   // conduction spin
   Sc("singlet", "triplet")   = -std::sqrt(0.75);
   Sc("triplet", "singlet")   =  0.5;
   Sc("triplet", "triplet")   =  std::sqrt(0.5);

   // local spin
   Sf("singlet", "triplet")   =  std::sqrt(0.75);
   Sf("triplet", "singlet")   = -0.5;
   Sf("triplet", "triplet")   = std::sqrt(0.5);
   Sf("holon",   "holon")     =  std::sqrt(0.75);

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
   Site[Sym1] = N;
   Site["Hu"] = Hu;
   Site[Sym2] = S;
   Site[Sym2+'f'] = Sf;
   Site[Sym2+'c'] = Sc;
   Site["C_A"] = C_A;
   Site["CH_A"] = CH_A;
   Site["CP_A"] = CP_A;
   Site["CHP_A"] = CHP_A;
   Site["C_B"] = C_B;
   Site["CH_B"] = CH_B;
   Site["CP_B"] = CP_B;
   Site["CHP_B"] = CHP_B;
   Site["N_S"] = Ns;
   Site["N_H"] = Nh;
   Site[Sym2+'c'+Sym2+'f'] = ScSf;
   return Site;
}

inline
LatticeSite CreateSO4KondoSiteA(std::string const& Sym1 = "Q", std::string const& Sym2 = "S")
{
   LatticeSite Site = CreateSO4KondoSiteCommon(Sym1, Sym2);
   Site["C"] = Site["C_A"];
   Site["CH"] = Site["CH_A"];
   Site["CP"] = Site["CP_A"];
   Site["CHP"] = Site["CHP_A"];
   return Site;
}

inline
LatticeSite CreateSO4KondoSiteB(std::string const& Sym1 = "Q", std::string const& Sym2 = "S")
{
   LatticeSite Site = CreateSO4KondoSiteCommon(Sym1, Sym2);
   Site["C"] = Site["C_B"];
   Site["CH"] = Site["CH_B"];
   Site["CP"] = Site["CP_B"];
   Site["CHP"] = Site["CHP_B"];
   return Site;
}
