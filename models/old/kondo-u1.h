// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/old/kondo-u1.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
LatticeSite CreateU1KondoSite(std::string const& Sym = "N")
{
   SymmetryList Symmetry(Sym+":U(1)");
   QuantumNumbers::QNConstructor<QuantumNumbers::U1> QN(Symmetry);
   SiteBasis Basis(Symmetry);
   SiteOperator Cup, Cdown, CHup, CHdown, P, R, N, Sp, Sm, Sz, Sfp, Sfm, Sfz, Scp, Scm, Scz,
      ScpSfm, ScmSfp, SczSfz, ScSf, I, Hu, Pdouble, 
      Ep, Em, Ns, Nh, Pg, CupP, CdownP, CHupP, CHdownP, Rs, Sx, Sy, Sfx, Sfy, Scx, Scy,
      mSz, Top;
   LatticeSite Site;

   Basis.push_back("empty-up",  QN(0));
   Basis.push_back("empty-down",  QN(0));
   Basis.push_back("double-up", QN(2));
   Basis.push_back("double-down", QN(2));
   Basis.push_back("singlet", QN(1));            // (1/sqrt(2)) ( |up DOWN> - |down UP>)
   Basis.push_back("triplet-down", QN(1));
   Basis.push_back("triplet-zero", QN(1));
   Basis.push_back("triplet-up", QN(1));

   Cup = SiteOperator(Basis, QN(-1), LatticeCommute::Fermionic);
   P = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   R = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   N = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   Hu = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   Pdouble = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   Pg = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   Sp = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   Sz = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   mSz = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   Top = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   Sfp = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   Sfz = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   Scp = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   Scz = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   I = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   Ns = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   Nh = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   Rs = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   ScSf = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);

   // annihilate up fermion
   Cup("empty-down",   "singlet")      =  sqrt(0.5);
   Cup("empty-down",   "triplet-zero") =  sqrt(0.5);
   Cup("empty-up",     "triplet-up")   =  1.0;
   Cup("singlet",      "double-up")    = -sqrt(0.5);
   Cup("triplet-zero", "double-up")    =  sqrt(0.5);
   Cup("triplet-down", "double-down")  =  1.0;

   // spin reflection
   Rs("empty-down",   "empty-up")     =  1.0;
   Rs("empty-up",     "empty-down")   =  1.0;

   Rs("double-down",  "double-up")    = -1.0;   // double-occupied singlet is antisymmetric
   Rs("double-up",    "double-down")  = -1.0;
   
   Rs("singlet",      "singlet")      = -1.0; // again, an antisymmetric triplet
   
   Rs("triplet-down", "triplet-up")   =  1.0;
   Rs("triplet-zero", "triplet-zero") =  1.0;
   Rs("triplet-up",   "triplet-down") =  1.0;

   // annihilate down fermion
   Cdown = Rs * Cup * Rs;

   // create up fermion
   CHup = adjoint(Cup);

   // create down fermion
   CHdown = adjoint(Cdown);

   TRACE(Cup)(CHup)(Cdown)(CHdown);

   // create double occupied site (eta pairing)
   Ep = CHup * CHdown;
   Em = Cdown * Cup;

   // parity = (-1)^N
   P("empty-up",     "empty-up")     =  1;
   P("empty-down",   "empty-down")   =  1;
   P("singlet",      "singlet")      = -1;
   P("triplet-down", "triplet-down") = -1;
   P("triplet-zero", "triplet-zero") = -1;
   P("triplet-up",   "triplet-up")   = -1;
   P("double-up",    "double-up")    =  1;
   P("double-down",  "double-down")  =  1;

   // parity modified creation/annihilation
   CupP = Cup * P;
   CdownP = Cdown * P;
   CHupP = CHup * P;
   CHdownP = CHdown * P;

   // spatial reflection
   R("empty-up",     "empty-up")     =  1;
   R("empty-down",   "empty-down")   =  1;
   R("singlet",      "singlet")      =  1;
   R("triplet-down", "triplet-down") =  1;
   R("triplet-zero", "triplet-zero") =  1;
   R("triplet-up",   "triplet-up")   =  1;
   R("double-up",    "double-up")    = -1;
   R("double-down",  "double-down")  = -1;

   // particle number
   N("singlet",      "singlet")      =  1;
   N("triplet-down", "triplet-down") =  1;
   N("triplet-zero", "triplet-zero") =  1;
   N("triplet-up",   "triplet-up")   =  1;
   N("double-up",    "double-up")    =  2;
   N("double-down",  "double-down")  =  2;

   // symmetrized coulomb operator = (n_up - 1/2) * (n_down - 1/2)
   Hu("empty-up",     "empty-up")     =  0.25;
   Hu("empty-down",   "empty-down")   =  0.25;
   Hu("singlet",      "singlet")      = -0.25;
   Hu("triplet-down", "triplet-down") = -0.25;
   Hu("triplet-zero", "triplet-zero") = -0.25;
   Hu("triplet-up",   "triplet-up")   = -0.25;
   Hu("double-up",    "double-up")    =  0.25;
   Hu("double-down",  "double-down")  =  0.25;

   // projection onto double occupied state == asymmetric coulomb operator == n_up * n_down
   Pdouble("double-up",   "double-up")   = 1.0;
   Pdouble("double-down", "double-down") = 1.0;

   // identity
   I("empty-up",     "empty-up")     =  1;
   I("empty-down",   "empty-down")   =  1;
   I("singlet",      "singlet")      =  1;
   I("triplet-down", "triplet-down") =  1;
   I("triplet-zero", "triplet-zero") =  1;
   I("triplet-up",   "triplet-up")   =  1;
   I("double-up",    "double-up")    =  1;
   I("double-down",  "double-down")  =  1;

   Pg = I - Pdouble; // Gutzwiller projector
   
   // local spin
   Sfz("empty-up",     "empty-up")     =  0.5;
   Sfz("empty-down",   "empty-down")   = -0.5;
   Sfz("triplet-zero", "singlet")      = -0.5;
   Sfz("singlet",      "triplet-zero") = -0.5;
   Sfz("triplet-up",   "triplet-up")   =  0.5;
   Sfz("triplet-down", "triplet-down") = -0.5;
   Sfz("empty-up",     "empty-up")     =  0.5;
   Sfz("empty-down",   "empty-down")   = -0.5;

   Sfp("empty-up",     "empty-down")   = 1.0;
   Sfp("double-up",    "double-down")  = 1.0;
   Sfp("singlet",      "triplet-down") = -sqrt(0.5);
   Sfp("triplet-zero", "triplet-down") =  sqrt(0.5);
   Sfp("triplet-up",   "singlet")      =  sqrt(0.5);
   Sfp("triplet-up",   "triplet-zero") =  sqrt(0.5);

   Sfm = adjoint(Sfp);

   //   CHECK_CLOSE(norm_frob(Sfm - Rs*Sfp*Rs), 0.0)(Sfm)(Rs*Sfp*Rs);

   // conduction spin
   Scz("triplet-zero", "singlet")      =  0.5;
   Scz("singlet",      "triplet-zero") =  0.5;
   Scz("triplet-up",   "triplet-up")   =  0.5;
   Scz("triplet-down", "triplet-down") = -0.5;

   Scp("singlet",      "triplet-down") =  sqrt(0.5);
   Scp("triplet-zero", "triplet-down") =  sqrt(0.5);
   Scp("triplet-up",   "singlet")      = -sqrt(0.5);
   Scp("triplet-up",   "triplet-zero") =  sqrt(0.5);

   Scm = adjoint(Scp);

   //   CHECK_CLOSE(norm_frob(Scm - Rs*Scp*Rs), 0.0)(Scm)(Rs*Scp*Rs);

   // total spin
   Sz("empty-up",     "empty-up")     =  0.5;
   Sz("empty-down",   "empty-down")   = -0.5;
   Sz("triplet-down", "triplet-down") = -1.0;
   Sz("triplet-up",   "triplet-up")   =  1.0;
   Sz("double-up",    "double-up")    =  0.5;
   Sz("double-down",  "double-down")  = -0.5;

   Sp("empty-up",      "empty-down")  =  1.0;
   Sp("double-up",     "double-down") =  1.0;
   Sp("triplet-zero", "triplet-down") =  sqrt(2.0);
   Sp("triplet-up",   "triplet-zero") =  sqrt(2.0);

   Sm = adjoint(Sp);

   // (-1)^Sz - string order parameter
   mSz("empty-up",     "empty-up")     =  std::complex<double>(0.0,  1.0);
   mSz("empty-down",   "empty-down")   =  std::complex<double>(0.0, -1.0);
   mSz("singlet",      "singlet")      =  1.0;
   mSz("triplet-zero", "triplet-zero") =  1.0;
   mSz("triplet-down", "triplet-down") = -1.0;
   mSz("triplet-up",   "triplet-up")   = -1.0;
   mSz("double-up",    "double-up")    =  std::complex<double>(0.0,  1.0);
   mSz("double-down",  "double-down")  =  std::complex<double>(0.0, -1.0);

   // Modified string operator for the doped topological insulator
   Top("empty-up",     "empty-up")     =  std::complex<double>(0.0,  1.0);
   Top("empty-down",   "empty-down")   =  std::complex<double>(0.0, -1.0);
   Top("singlet",      "singlet")      = -1.0;
   Top("triplet-zero", "triplet-zero") =  1.0;
   Top("triplet-down", "triplet-down") = -1.0;
   Top("triplet-up",   "triplet-up")   = -1.0;
   Top("double-up",    "double-up")    =  std::complex<double>(0.0,  1.0);
   Top("double-down",  "double-down")  =  std::complex<double>(0.0, -1.0);


   //   CHECK_CLOSE(norm_frob(Sp - (Scp + Sfp)), 0.0)(Sp)(Scp+Sfp);

   // number of spins in the conduction band
   Ns("singlet",      "singlet")      = 1;
   Ns("triplet-down", "triplet-down") = 1;
   Ns("triplet-zero", "triplet-zero") = 1;
   Ns("triplet-up",   "triplet-up")   = 1;

   // Local spin-spin interaction
   ScpSfm = Scp * Sfm;
   ScmSfp = Scm * Sfp;
   SczSfz = Scz * Sfz;

   // Isotropic spin-spin interaction
   ScSf("singlet",      "singlet")      = -0.75;
   ScSf("triplet-down", "triplet-down") =  0.25;
   ScSf("triplet-zero", "triplet-zero") =  0.25;
   ScSf("triplet-up",   "triplet-up")   =  0.25;

   // the x,y components of the spin
   // Sp = Sx + iSy
   // Sm = Sx - iSy
   // so that 
   // Sx = (Sp + Sm) / 2
   // Sy = (Sp - Sm)/2i

   std::complex<double> const i(0.0, 1.0);

   Sx = 0.5*(Sp + Sm);
   Sy = (-0.5*i)*(Sp - Sm);

   Scx = 0.5*(Scp + Scm);
   Scy = (-0.5*i)*(Scp - Scm);

   Sfx = 0.5*(Sfp + Sfm);
   Sfy = (-0.5*i)*(Sfp - Sfm);

   // number of holons
   Nh = I - Ns;

   Site["I"] = I;
   Site["P"] = P;
   Site[Sym] = N;
   Site["Hu"] = Hu;
   Site["Pdouble"] = Pdouble;
   Site["Pg"] = Pg;
   Site["Sp"] = Sp;
   Site["Sm"] = Sm;
   Site["Sz"] = Sz;
   Site["mSz"] = mSz;
   Site["Top"] = Top;
   Site["Sfp"] = Sfp;
   Site["Sfm"] = Sfm;
   Site["Sfz"] = Sfz;
   Site["Scp"] = Scp;
   Site["Scm"] = Scm;
   Site["Scz"] = Scz;
   Site["Cup"] = Cup;
   Site["CHup"] = CHup;
   Site["Cdown"] = Cdown;
   Site["CHdown"] = CHdown;
   Site["CupP"] = CupP;
   Site["CHupP"] = CHupP;
   Site["CdownP"] = CdownP;
   Site["CHdownP"] = CHdownP;
   Site["Ep"] = Ep;
   Site["Em"] = Em;
   Site["N_S"] = Ns;
   Site["N_H"] = Nh;
   Site["ScpSfm"] = ScpSfm;
   Site["ScmSfp"] = ScmSfp;
   Site["SczSfz"] = SczSfz;
   Site["ScSf"] = ScSf;
   Site["Sx"] = Sx;
   Site["Sy"] = Sy;
   Site["Sfx"] = Sfx;
   Site["Sfy"] = Sfy;
   Site["Scx"] = Scx;
   Site["Scy"] = Scy;

   return Site;
}
