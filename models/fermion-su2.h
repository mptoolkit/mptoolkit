// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/fermion-su2.h
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
LatticeSite FermionSU2(std::string const& Sym2 = "S")
{
   SymmetryList Symmetry(Sym2+":SU(2)");
   QuantumNumbers::QNConstructor<QuantumNumbers::SU2> QN(Symmetry);
   SiteBasis Basis(Symmetry);
   SiteOperator C, CH, P, R, N, S, I, Hu, Pdouble, Qp, Qm, Qz, N_S, N_H, Pg, ES;
   LatticeSite Site("SU(2) Fermion");

   Basis.push_back("empty",  QN(0));
   Basis.push_back("double", QN(0));
   Basis.push_back("single", QN(0.5));

   OperatorDescriptions OpDescriptions;
   OpDescriptions.add_operators()
      ("I"       , "identity")
      ("R"       , "reflection")
      ("P"       , "fermion parity")
      ("S"       , "spin vector operator")
      ("C"       , "annihilation operator")
      ("CH"      , "creation operator")
      ("N"       , "number operator")
      ("N_S"     , "number of spins")
      ("N_H"     , "number of holons")
      ("Qp"      , "eta raising operator (create double-occupied site)")
      ("Qm"      , "eta lowering operator (annihiliate double-occupied site)")
      ("Qz"      , "eta z operator, equivalent to (N-1)/2")
      ("ES"      , "exp(i*pi*s)")
      ("Hu"      , "symmetrized Coulomb operator (n_up - 1/2) * (n_down - 1/2)")
      ("Pdouble" , "projector onto the double-occupied site") 
      ("Pg"      , "Gutswiller projector = 1-Pdouble")
      ;

   C = SiteOperator(Basis, QN(0.5), LatticeCommute::Fermionic);
   CH = SiteOperator(Basis, QN(0.5), LatticeCommute::Fermionic);
   P = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   R = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   N = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   Hu = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   Pdouble = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   Pg = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   S = SiteOperator(Basis, QN(1), LatticeCommute::Bosonic);
   I = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   N_S = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   N_H = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);

   // annihilate fermion
   C("empty",  "single")    =  std::sqrt(2.0);
   C("single", "double")    =  1;
   
   // create fermion
   CH = adjoint(C);

   // create double occupied site (eta pairing)
   Qp = sqrt(2.0) * prod(CH, CH, QN(0));
   Qm = sqrt(2.0) * prod(C, C, QN(0));

   // parity = (-1)^N   
   P("empty",  "empty")     =  1;
   P("single", "single")    = -1;
   P("double", "double")    =  1;

   // spatial reflection   
   R("empty",  "empty")     =  1;
   R("single", "single")    =  1;
   R("double", "double")    = -1;
 
   // particle number  
   N("single", "single")    =  1;
   N("double", "double")    =  2;

   // symmetrized coulomb operator = (n_up - 1/2) * (n_down - 1/2)
   Hu("empty",  "empty")     =  0.25;
   Hu("single", "single")    = -0.25;
   Hu("double", "double")    =  0.25;

   // projection onto double occupied state == asymmetric coulomb operator == n_up * n_down
   Pdouble("double", "double") = 1.0;

   // identity
   I("empty",  "empty")     =  1;
   I("single", "single")    =  1;
   I("double", "double")    =  1;

   Pg = I - Pdouble; // Gutzwiller projector
   
   // S
   S("single", "single")   = std::sqrt(0.75);

   // exp(i * pi * S)
   ES = I;
   ES("single", "single") = std::complex<double>(0.0, 1.0);

   // number of spins
   N_S("single", "single") = 1;

   // number of holons
   N_H = I - N_S;

   Site["I"] = I;
   Site["P"] = P;
   Site["N"] = N;
   Site["Hu"] = Hu;
   Site["Pdouble"] = Pdouble;
   Site["Pg"] = Pg;
   Site[Sym2] = S;
   Site["C"] = C;
   Site["CH"] = CH;
   Site["Qp"] = Qp;
   Site["Qm"] = Qm;
   Site["Qz"] = 0.5*(N-I);
   Site["N_S"] = N_S;
   Site["N_H"] = N_H;
   Site["ES"] = ES;
   Site["R"] = R;

   Site.set_operator_descriptions(OpDescriptions);

   return Site;
}
