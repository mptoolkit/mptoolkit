// -*- C++ -*- $Id: hubbard-u1u1.h 1350 2014-04-01 17:26:24Z ianmcc $

#include "lattice/latticesite.h"
#include "quantumnumbers/null-quantumnumber.h"

inline
LatticeSite CreateHubbardSite()
{
   using QuantumNumbers::Null;
   SiteBasis Basis(QuantumNumbers::NullSymmetryList);
   SiteOperator CHup, Cup, CHdown, Cdown, P, R, N, Sp, Sm, Sz, Qp, Qm, Qz, I, Hu, Pdouble, X, N_S, N_H, ES, EQ,
      Sx, Sy, mSx, mSy;
   LatticeSite Site;

   Basis.push_back("empty",  Null);
   Basis.push_back("double", Null);
   Basis.push_back("down" ,  Null);
   Basis.push_back("up" ,    Null);

   CHup = SiteOperator(Basis, Null, LatticeCommute::Fermionic);
   Cup = SiteOperator(Basis, Null, LatticeCommute::Fermionic);
   CHdown = SiteOperator(Basis, Null, LatticeCommute::Fermionic);
   Cdown = SiteOperator(Basis, Null, LatticeCommute::Fermionic);
   P = SiteOperator(Basis, Null, LatticeCommute::Bosonic);
   R = SiteOperator(Basis, Null, LatticeCommute::Bosonic);
   N = SiteOperator(Basis, Null, LatticeCommute::Bosonic);
   Pdouble = SiteOperator(Basis, Null, LatticeCommute::Bosonic);
   Hu = SiteOperator(Basis, Null, LatticeCommute::Bosonic);
   Sp = SiteOperator(Basis, Null, LatticeCommute::Bosonic);
   Sm = SiteOperator(Basis, Null, LatticeCommute::Bosonic);
   Sz = SiteOperator(Basis, Null, LatticeCommute::Bosonic);
   Sx = SiteOperator(Basis, Null, LatticeCommute::Bosonic);
   Sy = SiteOperator(Basis, Null, LatticeCommute::Bosonic);
   Qp = SiteOperator(Basis, Null, LatticeCommute::Bosonic);
   Qm = SiteOperator(Basis, Null, LatticeCommute::Bosonic);
   Qz = SiteOperator(Basis, Null, LatticeCommute::Bosonic);
   I = SiteOperator(Basis, Null, LatticeCommute::Bosonic);
   X = SiteOperator(Basis, Null, LatticeCommute::Bosonic);
   N_S = SiteOperator(Basis, Null, LatticeCommute::Bosonic);
   N_H = SiteOperator(Basis, Null, LatticeCommute::Bosonic);

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

   X("empty",  "empty")     = 1;
   X("up",     "up")        = std::complex<double>(0,1);
   X("down",   "down")      = std::complex<double>(0,1);
   X("double", "double")    = 1;

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

   Sx = Sp + Sm;
   Sy = std::complex<double>(0.0,  1.0) * (Sp - Sm);

   mSx =  LinearAlgebra::exp(std::complex<double>(0.0,math_const::pi) * Sx);
   mSy =  LinearAlgebra::exp(std::complex<double>(0.0,math_const::pi) * Sy);

   // Q^+
   Qp("double", "empty") = 1;

   // Q^-
   Qm("empty", "double") = 1;

   // Q^z
   Qz("double", "double") =  0.5;
   Qz("empty",  "empty")  = -0.5;
   
   // exp(i * pi * Sz)
   ES = I;
   ES("up",   "up")   = std::complex<double>(0.0,  1.0);
   ES("down", "down") = std::complex<double>(0.0, -1.0);

   // exp(i * pi * Qz)
   EQ = I;
   EQ("double", "double") = std::complex<double>(0.0,  1.0);
   EQ("empty",  "empty")  = std::complex<double>(0.0, -1.0);

   // Same as ES above
   Site["mSz"] = LinearAlgebra::exp(std::complex<double>(0.0,math_const::pi) * Sz);
   Site["m2Sz"] = LinearAlgebra::exp(std::complex<double>(0.0,2.0*math_const::pi) * Sz);
   Site["Qx"] = Qp + Qm;
   Site["Qy"] = std::complex<double>(0.0, 1.0) * (Qp - Qm);
   Site["I"] = I;
   Site["P"] = P;
   Site["N"] = N;
   Site["Hu"] = Hu;
   Site["Pdouble"] = Pdouble;
   Site["Sp"] = Sp;
   Site["Sm"] = Sm;
   Site["Sz"] = Sz;
   Site["Sx"] = Sx;
   Site["Sy"] = Sy;
   Site["mSx"] = mSx;
   Site["mSy"] = mSy;
   Site["Qp"] = Qp;
   Site["Qm"] = Qm;
   Site["Qz"] = Qz;
   Site["R"] = R;
   Site["CHup"] = CHup;
   Site["Cup"] = Cup;
   Site["CHdown"] = CHdown;
   Site["Cdown"] = Cdown;
   Site["CHupP"] = CHup*P;
   Site["CupP"] = Cup*P;
   Site["CHdownP"] = CHdown*P;
   Site["CdownP"] = Cdown*P;
   Site["X"] = X;
   Site["N_S"] = N_S;
   Site["N_H"] = N_H;
   Site["ES"] = ES;
   Site["EQ"] = EQ;
   return Site;
}
