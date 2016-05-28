// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/fermion-u1u1.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "lattice/latticesite.h"
#include "quantumnumbers/u1.h"

inline
LatticeSite FermionU1U1(std::string const& Sym1 = "N", 
			std::string const& Sym2 = "Sz",
			std::string const& ParityOp = "P")
{
   SymmetryList Symmetry(Sym1+":U(1),"+Sym2+":U(1)");
   QuantumNumbers::QNConstructor<QuantumNumbers::U1,QuantumNumbers::U1> QN(Symmetry);
   SiteBasis Basis(Symmetry);
   SiteOperator CHup, Cup, CHdown, Cdown, P, R, N, Sp, Sm, Sz, Qp, Qm, Qz, I, Nup, Ndown,
      Hu, Pdouble, ES, N_S, N_H;
   LatticeSite Site("U(1)xU(1) Fermion");

   Basis.push_back("empty",  QN(0, 0));
   Basis.push_back("double", QN(2, 0));
   Basis.push_back("down" ,  QN(1,-0.5));
   Basis.push_back("up" ,    QN(1, 0.5));

   OperatorDescriptions OpDescriptions;
   OpDescriptions.add_operators()
      ("I"       , "identity")
      ("R"       , "reflection")
      ("P"       , "fermion parity")
      ("Sp"      , "spin raising operator")
      ("Sm"      , "spin lowering operator")
      ("Sz"      , "z-component of spin")
      ("Cup"     , "annihilation up spin fermion")
      ("Cdown"   , "annihilation down spin fermion")
      ("CHup"    , "creation up spin fermion")
      ("CHdown"  , "creation down spin fermion")
      ("Qp"      , "eta raising operator (create double-occupied site)")
      ("Qm"      , "eta lowering operator (annihiliate double-occupied site)")
      ("Qz"      , "eta z-operator, equivalent to (N-1)/2")
      ("N"       , "number operator")
      ("Nup"     , "number of up spins")
      ("Ndown"   , "number of down spins")
      ("N_S"     , "number of spins")
      ("N_H"     , "number of holons")
      ("ES"      , "exp(i*pi*s)")
      ("Hu"      , "symmetrized Coulomb operator (n_up - 1/2) * (n_down - 1/2)")
      ("Pdouble" , "projector onto the double-occupied site") 
      ("Pg"      , "Gutswiller projector = 1-Pdouble")
      ;


   LatticeCommute Fermionic = ParityOp == "P" ? LatticeCommute::Fermionic : LatticeCommute(ParityOp);

   CHup = SiteOperator(Basis, QN(1,0.5), Fermionic);
   Cup = SiteOperator(Basis, QN(-1,-0.5), Fermionic);
   CHdown = SiteOperator(Basis, QN(1,-0.5), Fermionic);
   Cdown = SiteOperator(Basis, QN(-1,0.5), Fermionic);
   P = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   R = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   N = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   Nup = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   Ndown = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   Pdouble = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   Hu = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   Sp = SiteOperator(Basis, QN(0,1), LatticeCommute::Bosonic);
   Sm = SiteOperator(Basis, QN(0,-1), LatticeCommute::Bosonic);
   Sz = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   Qp = SiteOperator(Basis, QN(2,0), LatticeCommute::Bosonic);
   Qm = SiteOperator(Basis, QN(-2,0), LatticeCommute::Bosonic);
   Qz = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   I = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   ES = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   N_S = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   N_H = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);

   // create up spin
   CHup("up",     "empty")     =  1;
   CHup("double", "down")      =  1;
   
   // annihilate up spin
   Cup = adjoint(CHup);

   // create down spin   
   CHdown("down",   "empty")   =  1;
   CHdown("double", "up")      = -1;
   
   // annihilate down spin
   Cdown = adjoint(CHdown);

   // parity = (-1)^N   
   P("empty",  "empty")     =  1;
   P("up",     "up")        = -1;
   P("down",   "down")      = -1;
   P("double", "double")    =  1;

   ES("empty",  "empty")     = 1;
   ES("up",     "up")        = std::complex<double>(0,1);
   ES("down",   "down")      = std::complex<double>(0,1);
   ES("double", "double")    = 1;

   N_H("empty",  "empty")     = 1;
   N_H("double", "double")    = 1;

   N_S("up",     "up")        = 1;
   N_S("down",   "down")      = 1;

   // spatial reflection   
   R("empty",  "empty")     =  1;
   R("up",     "up")        =  1;
   R("down",   "down")      =  1;
   R("double", "double")    = -1;
 
   // particle number  
   N("up",     "up")        =  1;
   N("down",   "down")      =  1;
   N("double", "double")    =  2;

   Nup("up",     "up")       =  1;
   Nup("double", "double")   =  1;

   Ndown("down",   "down")   =  1;
   Ndown("double", "double") =  1;

   // symmetrized coulomb operator = (n_up - 1/2) * (n_down - 1/2)
   Hu("empty",  "empty")     =  0.25;
   Hu("up",     "up")        = -0.25;
   Hu("down",   "down")      = -0.25;
   Hu("double", "double")    =  0.25;

   // projection onto double occupied state 
   // == asymmetric coulomb operator == n_up * n_down
   Pdouble("double", "double") = 1.0;

   // identity
   I("empty",  "empty")     =  1;
   I("up",     "up")        =  1;
   I("down",   "down")      =  1;
   I("double", "double")    =  1;
   
   // S^+
   Sp("up",   "down")      = 1;

   // S^-
   Sm("down", "up")        = 1;
   
   // z-component of spin
   Sz("up",   "up")       =  0.5;
   Sz("down", "down")     = -0.5;

   // Q^+
   Qp("double", "empty") = 1;

   // Q^-
   Qm("empty", "double") = 1;

   // Q^z
   Qz("double", "double") =  0.5;
   Qz("empty",  "empty")  = -0.5;
   
   Site["I"] = I;
   Site[ParityOp] = P;
   Site["N"] = N;
   Site["Nup"] = Nup;
   Site["Ndown"] = Ndown;
   Site[Sym1] = N;
   Site["Hu"] = Hu;
   Site["Pdouble"] = Pdouble;
   Site["Pg"] = I - Pdouble;
   Site["Sp"] = Sp;
   Site["Sm"] = Sm;
   Site["Sz"] = Sz;
   Site[Sym2] = Sz;
   Site["Qp"] = Qp;
   Site["Qm"] = Qm;
   Site["Qz"] = Qz;
   Site["R"] = R;
   Site["CHup"] = CHup;
   Site["Cup"] = Cup;
   Site["CHdown"] = CHdown;
   Site["Cdown"] = Cdown;
   Site["ES"] = ES;
   Site["N_S"] = N_S;
   Site["N_H"] = N_H;

   Site.set_operator_descriptions(OpDescriptions);

   return Site;
}
