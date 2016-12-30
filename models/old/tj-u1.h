// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/old/tj-u1.h
//
// Copyright (C) 2015-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#if !defined(TJ_U1_H_SDHFOUIERYT8934758934789)
#define TJ_U1_H_SDHFOUIERYT8934758934789

#include "lattice/latticesite.h"
#include "quantumnumbers/u1.h"




inline
LatticeSite CreateU1tJSite(std::string const& Sym1 = "N")
{
   SymmetryList Symmetry(Sym1+":U(1)");
   QuantumNumbers::QNConstructor<QuantumNumbers::U1> QN(Symmetry);
   SiteBasis Basis(Symmetry);
   SiteOperator CHup, Cup, CHdown, Cdown, P, R, N, Sp, Sm, Sz, Qp, Qm, Qz, I, Hu, Pdouble, X, N_S, N_H;
   LatticeSite Site;

   Basis.push_back("empty",  QN(0));
   Basis.push_back("down" ,  QN(1));
   Basis.push_back("up" ,    QN(1));

   CHup = SiteOperator(Basis, QN(1), LatticeCommute::Fermionic);
   Cup = SiteOperator(Basis, QN(-1), LatticeCommute::Fermionic);
   CHdown = SiteOperator(Basis, QN(1), LatticeCommute::Fermionic);
   Cdown = SiteOperator(Basis, QN(-1), LatticeCommute::Fermionic);
   P = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   R = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   N = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   Hu = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   Sp = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   Sm = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   Sz = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   I = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   X = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   N_S = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   N_H = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);

   // create up spin
   CHup("up",     "empty")     =  1;

   // annihilate up spin
   Cup = adjoint(CHup);

   // create down spin
   CHdown("down",   "empty")   =  1;

   // annihilate down spin
   Cdown = adjoint(CHdown);

   // parity = (-1)^N
   P("empty",  "empty")     =  1;
   P("up",     "up")        = -1;
   P("down",   "down")      = -1;

   X("empty",  "empty")     = 1;
   X("up",     "up")        = std::complex<double>(0,1);
   X("down",   "down")      = std::complex<double>(0,1);

   N_H("empty",  "empty")     = 1;

   N_S("up",     "up")        = 1;
   N_S("down",   "down")      = 1;

   // spatial reflection
   R("empty",  "empty")     =  1;
   R("up",     "up")        =  1;
   R("down",   "down")      =  1;

   // particle number
   N("up",     "up")        =  1;
   N("down",   "down")      =  1;

   // symmetrized coulomb operator = (n_up - 1/2) * (n_down - 1/2)
   Hu("empty",  "empty")     =  0.25;
   Hu("up",     "up")        = -0.25;
   Hu("down",   "down")      = -0.25;

   // identity
   I("empty",  "empty")     =  1;
   I("up",     "up")        =  1;
   I("down",   "down")      =  1;

   // S^+
   Sp("up",   "down")      = 1;

   // S^-
   Sm("down", "up")        = 1;

   // z-component of spin
   Sz("up",   "up")       =  0.5;
   Sz("down", "down")     = -0.5;

   Site["I"] = I;
   Site["P"] = P;
   Site["N"] = N;
   Site[Sym1] = N;
   Site["Hu"] = Hu;
   Site["Pdouble"] = Pdouble;
   Site["Sp"] = Sp;
   Site["Sm"] = Sm;
   Site["Sz"] = Sz;
   Site["R"] = R;
   Site["CHup"] = CHup;
   Site["Cup"] = Cup;
   Site["CHdown"] = CHdown;
   Site["Cdown"] = Cdown;
   Site["CHupP"] = prod(CHup, P, QN(1));
   Site["CupP"] = prod(Cup, P, QN(-1));
   Site["CHdownP"] = prod(CHdown, P, QN(1));
   Site["CdownP"] = prod(Cdown, P, QN(-1));
   Site["X"] = X;
   Site["N_S"] = N_S;
   Site["N_H"] = N_H;
   return Site;
}

#endif
