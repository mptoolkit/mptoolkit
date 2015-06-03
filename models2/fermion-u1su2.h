// -*- C++ -*- $Id$

#include "lattice/latticesite.h"
#include "quantumnumbers/u1.h"
#include "quantumnumbers/su2.h"

inline
LatticeSite FermionU1SU2(std::string const& Sym1 = "N", std::string const& Sym2 = "S")
{
   SymmetryList Symmetry(Sym1+":U(1),"+Sym2+":SU(2)");
   QuantumNumbers::QNConstructor<QuantumNumbers::U1,QuantumNumbers::SU2> QN(Symmetry);
   SiteBasis Basis(Symmetry);
   SiteOperator C, CH, P, R, N, S, I, Hu, Pdouble, Ep, Em, Ns, Nh, Pg, CP, CHP, ES;
   LatticeSite Site;

   Basis.push_back("empty",  QN(0, 0));
   Basis.push_back("double", QN(2, 0));
   Basis.push_back("single", QN(1, 0.5));

   C = SiteOperator(Basis, QN(-1,0.5), LatticeCommute::Fermionic);
   CH = SiteOperator(Basis, QN(1,0.5), LatticeCommute::Fermionic);
   P = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   R = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   N = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   Hu = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   Pdouble = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   Pg = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   S = SiteOperator(Basis, QN(0,1), LatticeCommute::Bosonic);
   I = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   Ns = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   Nh = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);

   // annihilate fermion
   C("empty",  "single")    =  std::sqrt(2.0);
   C("single", "double")    =  1;
   
   // create fermion
   CH = adjoint(C);

   // create double occupied site (eta pairing)
   Ep = sqrt(2.0) * prod(CH, CH, QN(2,0));
   Em = sqrt(2.0) * prod(C, C, QN(-2,0));

   // parity = (-1)^N   
   P("empty",  "empty")     =  1;
   P("single", "single")    = -1;
   P("double", "double")    =  1;

   // parity modified creation/annihilation
   CP = prod(C, P, QN(-1,0.5));
   CHP = prod(CH, P, QN(1,0.5));

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
   Ns("single", "single") = 1;

   // number of holons
   Nh = I - Ns;

   Site["I"] = I;
   Site["P"] = P;
   Site[Sym1] = N;
   Site["Hu"] = Hu;
   Site["Pdouble"] = Pdouble;
   Site["Pg"] = Pg;
   Site[Sym2] = S;
   Site["C"] = C;
   Site["CH"] = CH;
   Site["CP"] = CP;
   Site["CHP"] = CHP;
   Site["Ep"] = Ep;
   Site["Em"] = Em;
   Site["N_S"] = Ns;
   Site["N_H"] = Nh;
   Site["ES"] = ES;
   Site["R"] = R;
   return Site;
}
