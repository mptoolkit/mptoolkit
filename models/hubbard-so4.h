// -*- C++ -*- $Id$

#if !defined(HUBBARD_SO4_H)
#define HUBBARD_SO4_H

#include "siteoperator/siteoperator.h"
#include "quantumnumbers/su2.h"
#include "siteoperator/block.h"
#include <cmath>

typedef Block<SiteOperator> SiteBlock;

using std::sqrt;

inline
SiteBlock CreateSO4HubbardSiteCommon(std::string const& Sym1 = "Q", std::string const& Sym2 = "S")
{
   SymmetryList Symmetry(Sym1+":SU(2),"+Sym2+":SU(2)");
   QuantumNumbers::QNConstructor<QuantumNumbers::SU2,QuantumNumbers::SU2> QN(Symmetry);
   SiteBasis Basis(Symmetry);
   SiteOperator C_A, C_B, CH_A, CH_B, P, R, N, N_S, N_H, S, Q, I, X;
   SiteBlock Block;

   Basis.push_back("holon",  QN(0.5, 0));
   Basis.push_back("spinon", QN(0, 0.5));

   C_A = SiteOperator(Basis, QN(0.5,0.5), LatticeCommute::Fermionic);
   C_B = SiteOperator(Basis, QN(0.5,0.5), LatticeCommute::Fermionic);
   P = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   R = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   N = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   N_S = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   N_H = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   S = SiteOperator(Basis, QN(0,1), LatticeCommute::Bosonic);
   Q = SiteOperator(Basis, QN(1,0), LatticeCommute::Bosonic);
   I = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   X = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);

   C_A("spinon", "holon")  = sqrt(2.0);
   C_A("holon",  "spinon") = sqrt(2.0);

   C_B("spinon", "holon")  =  sqrt(2.0);
   C_B("holon",  "spinon") = -sqrt(2.0);

   CH_A = adjoint(C_A);
   CH_B = adjoint(C_B);

   P("spinon", "spinon") = -1;
   P("holon",  "holon")  = 1;

   SiteOperator CP_A = prod(C_A, P, QN(0.5,0.5));
   SiteOperator CP_B = prod(C_B, P, QN(0.5,0.5));
   SiteOperator CHP_A = prod(CH_A, P, QN(0.5,0.5));
   SiteOperator CHP_B = prod(CH_B, P, QN(0.5,0.5));

   X("spinon", "spinon") = std::complex<double>(0, 1);
   X("holon", "holon") = 1;

   R("spinon", "spinon") = 1;
   R("holon",  "holon")  = 1;

   N_S("spinon", "spinon") = 1;

   N_H("holon", "holon") = 1;

   // this choice of signs makes the 0-component from the Wigner Eckart theorem 
   // equal to the usual definition of Sz.  Note that the "hermitian conjugate" of
   // S is S^\dagger = -S.  This arises from the signs picked up in taking the
   // transpose of a spin 1 operator.
   S("spinon", "spinon") = sqrt(0.75);

   Q("holon", "holon")   = sqrt(0.75);

   // identity
   I("spinon", "spinon") =  1;
   I("holon",  "holon")  =  1;
   
   Block["I"] = I;
   Block["P"] = P;
   Block["N"] = N;
   Block["N_S"] = N_S;
   Block["N_H"] = N_H;
   Block[Sym1] = Q;
   Block[Sym2] = S;
   Block["C_A"] = C_A;
   Block["C_B"] = C_B;
   Block["CH_A"] = CH_A;
   Block["CH_B"] = CH_B;
   Block["CP_A"] = CP_A;
   Block["CP_B"] = CP_B;
   Block["CHP_A"] = CHP_A;
   Block["CHP_B"] = CHP_B;
   Block["X"] = X;
   return Block;
}

inline
SiteBlock CreateSO4HubbardSiteA(std::string const& Sym1 = "Q", std::string const& Sym2 = "S")
{
   SiteBlock Site = CreateSO4HubbardSiteCommon(Sym1, Sym2);
   Site["C"] = Site["C_A"];
   Site["CH"] = Site["CH_A"];
   Site["CP"] = Site["CP_A"];
   Site["CHP"] = Site["CHP_A"];
   return Site;
}

inline
SiteBlock CreateSO4HubbardSiteB(std::string const& Sym1 = "Q", std::string const& Sym2 = "S")
{
   SiteBlock Site = CreateSO4HubbardSiteCommon(Sym1, Sym2);
   Site["C"] = Site["C_B"];
   Site["CH"] = Site["CH_B"];
   Site["CP"] = Site["CP_B"];
   Site["CHP"] = Site["CHP_B"];
   return Site;
}

#endif
