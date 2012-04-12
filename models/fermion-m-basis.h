// -*- C++ -*- $Id$
//
// A fermion site, designed for the nuclear shell model in the m-basis
// each site is hubbard-like, with empty, up,down and double-occupied states
// but the z-component of spin is a parameter.
//

#include "siteoperator/siteoperator.h"
#include "quantumnumbers/u1.h"
#include "siteoperator/block.h"

typedef Block<SiteOperator> SiteBlock;

inline
SiteBlock CreateFermionSite(half_int m, std::string const& Sym1 = "N", std::string const& Sym2 = "Sz")
{
   SymmetryList Symmetry(Sym1+":U(1),"+Sym2+":U(1)");
   QuantumNumbers::QNConstructor<QuantumNumbers::U1,QuantumNumbers::U1> QN(Symmetry);
   SiteBasis Basis(Symmetry);
   SiteOperator CHup, Cup, CHdown, Cdown, P, R, N, Sz, I, Pdouble, X, N_S, N_H;
   SiteBlock Block;

   Basis.push_back("empty",  QN(0, 0));
   Basis.push_back("up" ,    QN(1, m));
   Basis.push_back("down" ,  QN(1,-m));
   Basis.push_back("double", QN(2, 0));

   CHup = SiteOperator(Basis, QN(1,m), LatticeCommute::Fermionic);
   Cup = SiteOperator(Basis, QN(-1,-m), LatticeCommute::Fermionic);
   CHdown = SiteOperator(Basis, QN(1,-m), LatticeCommute::Fermionic);
   Cdown = SiteOperator(Basis, QN(-1,m), LatticeCommute::Fermionic);
   P = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   R = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   N = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   Pdouble = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   Sz = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
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

   // projection onto double occupied state 
   // == asymmetric coulomb operator == n_up * n_down
   Pdouble("double", "double") = 1.0;

   // identity
   I("empty",  "empty")     =  1;
   I("up",     "up")        =  1;
   I("down",   "down")      =  1;
   I("double", "double")    =  1;
   
   // z-component of spin
   Sz("up",   "up")       =  m.to_double();
   Sz("down", "down")     = -m.to_double();

   Block["I"] = I;
   Block["P"] = P;
   Block["N"] = N;
   Block[Sym1] = N;
   Block["Pdouble"] = Pdouble;
   Block["Sz"] = Sz;
   Block[Sym2] = Sz;
   Block["CHup"] = CHup;
   Block["Cup"] = Cup;
   Block["CHdown"] = CHdown;
   Block["Cdown"] = Cdown;
   Block["X"] = X;
   Block["N_S"] = N_S;
   Block["N_H"] = N_H;
   return Block;
}
