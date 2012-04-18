// -*- C++ -*- $Id$

#include "siteoperator/siteoperator.h"
#include "quantumnumbers/u1.h"
#include "quantumnumbers/su2.h"
#include "siteoperator/block.h"

typedef Block<SiteOperator> SiteBlock;

inline
SiteBlock CreateU1SU2tJSite(std::string const& Sym1 = "N", std::string const& Sym2 = "S")
{
   SymmetryList Symmetry(Sym1+":U(1),"+Sym2+":SU(2)");
   QuantumNumbers::QNConstructor<QuantumNumbers::U1,QuantumNumbers::SU2> QN(Symmetry);
   SiteBasis Basis(Symmetry);
   SiteOperator C, CH, P, R, N, S, I, Hu, Ns, Nh, Pg, CP, CHP;
   SiteBlock Block;

   Basis.push_back("empty",  QN(0, 0));
   Basis.push_back("single", QN(1, 0.5));

   C = SiteOperator(Basis, QN(-1,0.5), LatticeCommute::Fermionic);
   CH = SiteOperator(Basis, QN(1,0.5), LatticeCommute::Fermionic);
   P = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   R = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   N = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   Hu = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   S = SiteOperator(Basis, QN(0,1), LatticeCommute::Bosonic);
   I = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   Ns = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   Nh = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);

   // annihilate fermion
   C("empty",  "single")    =  std::sqrt(2.0);
   
   // create fermion
   CH = adjoint(C);

   // parity = (-1)^N   
   P("empty",  "empty")     =  1;
   P("single", "single")    = -1;

   // parity modified creation/annihilation
   CP = prod(C, P, QN(-1,0.5));
   CHP = prod(CH, P, QN(1,0.5));

   // spatial reflection   
   R("empty",  "empty")     =  1;
   R("single", "single")    =  1;
 
   // particle number  
   N("single", "single")    =  1;

   // symmetrized coulomb operator = (n_up - 1/2) * (n_down - 1/2)
   Hu("empty",  "empty")     =  0.25;
   Hu("single", "single")    = -0.25;

   // identity
   I("empty",  "empty")     =  1;
   I("single", "single")    =  1;

   // S
   S("single", "single")   = std::sqrt(0.75);

   // number of spins
   Ns("single", "single") = 1;

   // number of holons
   Nh = I - Ns;

   Block["I"] = I;
   Block["P"] = P;
   Block[Sym1] = N;
   Block["Hu"] = Hu;
   Block[Sym2] = S;
   Block["C"] = C;
   Block["CH"] = CH;
   Block["CP"] = CP;
   Block["CHP"] = CHP;
   Block["N_S"] = Ns;
   Block["N_H"] = Nh;
   return Block;
}
