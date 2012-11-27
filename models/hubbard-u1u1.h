// -*- C++ -*- $Id$

#include "siteoperator/siteoperator.h"
#include "quantumnumbers/u1.h"
#include "siteoperator/block.h"

typedef Block<SiteOperator> SiteBlock;

inline
SiteBlock CreateU1HubbardSite(std::string const& Sym1 = "N", 
			      std::string const& Sym2 = "Sz",
			      std::string const& ParityOp = "P")
{
   SymmetryList Symmetry(Sym1+":U(1),"+Sym2+":U(1)");
   QuantumNumbers::QNConstructor<QuantumNumbers::U1,QuantumNumbers::U1> QN(Symmetry);
   SiteBasis Basis(Symmetry);
   SiteOperator CHup, Cup, CHdown, Cdown, P, R, N, Sp, Sm, Sz, Qp, Qm, Qz, I, Hu, Pdouble, X, N_S, N_H, ES, EQ;
   SiteBlock Block;

#if defined(NEW_SITE_ORDERING)
#pragma warning "using new site ordering."
   Basis.push_back("empty",  QN(0, 0));
   Basis.push_back("double", QN(2, 0));
   Basis.push_back("down" ,  QN(1,-0.5));
   Basis.push_back("up" ,    QN(1, 0.5));
#elif defined(OLD_SITE_ORDERING)
   Basis.push_back("empty",  QN(0, 0));
   Basis.push_back("up" ,    QN(1, 0.5));
   Basis.push_back("down" ,  QN(1,-0.5));
   Basis.push_back("double", QN(2, 0));
#else
#pragma warning "using old site ordering."
   Basis.push_back("empty",  QN(0, 0));
   Basis.push_back("down" ,  QN(1,-0.5));
   Basis.push_back("up" ,    QN(1, 0.5));
   Basis.push_back("double", QN(2, 0));
#endif

   LatticeCommute Fermionic = ParityOp == "P" ? LatticeCommute::Fermionic : LatticeCommute(ParityOp);

   CHup = SiteOperator(Basis, QN(1,0.5), Fermionic);
   Cup = SiteOperator(Basis, QN(-1,-0.5), Fermionic);
   CHdown = SiteOperator(Basis, QN(1,-0.5), Fermionic);
   Cdown = SiteOperator(Basis, QN(-1,0.5), Fermionic);
   P = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   R = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   N = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   Pdouble = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   Hu = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   Sp = SiteOperator(Basis, QN(0,1), LatticeCommute::Bosonic);
   Sm = SiteOperator(Basis, QN(0,-1), LatticeCommute::Bosonic);
   Sz = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   Qp = SiteOperator(Basis, QN(2,0), LatticeCommute::Bosonic);
   Qm = SiteOperator(Basis, QN(-2,0), LatticeCommute::Bosonic);
   Qz = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   I = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   X = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
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

   Block["I"] = I;
   Block[ParityOp] = P;
   Block["N"] = N;
   Block[Sym1] = N;
   Block["Hu"] = Hu;
   Block["Pdouble"] = Pdouble;
   Block["Sp"] = Sp;
   Block["Sm"] = Sm;
   Block["Sz"] = Sz;
   Block[Sym2] = Sz;
   Block["Qp"] = Qp;
   Block["Qm"] = Qm;
   Block["Qz"] = Qz;
   Block["R"] = R;
   Block["CHup"] = CHup;
   Block["Cup"] = Cup;
   Block["CHdown"] = CHdown;
   Block["Cdown"] = Cdown;
   Block["CHupP"] = CHup*P;
   Block["CupP"] = Cup*P;
   Block["CHdownP"] = CHdown*P;
   Block["CdownP"] = Cdown*P;
   Block["X"] = X;
   Block["N_S"] = N_S;
   Block["N_H"] = N_H;
   Block["ES"] = ES;
   Block["EQ"] = EQ;
   return Block;
}
