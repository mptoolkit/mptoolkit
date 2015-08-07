// -*- C++ -*- $Id: fermion-so4.h 1495 2015-05-20 16:11:06Z ianmcc $

#if !defined(HUBBARD_SO4_H)
#define HUBBARD_SO4_H

#include "lattice/latticesite.h"
#include "quantumnumbers/su2.h"
#include <cmath>

inline
LatticeSite FermionSO4_Common(std::string const& Sym1 = "Q", std::string const& Sym2 = "S")
{
   SymmetryList Symmetry(Sym1+":SU(2),"+Sym2+":SU(2)");
   QuantumNumbers::QNConstructor<QuantumNumbers::SU2,QuantumNumbers::SU2> QN(Symmetry);
   SiteBasis Basis(Symmetry);
   SiteOperator C_A, C_B, CH_A, CH_B, P, R, N_S, N_H, S, Q, I, ES;
   LatticeSite Site("SO(4) Fermion");

   Basis.push_back("holon",  QN(0.5, 0));
   Basis.push_back("spinon", QN(0, 0.5));

   OperatorDescriptions OpDescriptions;
   OpDescriptions.add_operators()
      ("I"       , "identity")
      ("P"       , "fermion parity")
      (Sym2      , "spin vector operator")
      (Sym1      , "pseudo-spin vector operator")
      ("C_A"     , "annihilation operator for A site")
      ("CH_A"    , "creation operator for A site")
      ("C_B"     , "annihilation operator for B site")
      ("CH_B"    , "creation operator for B site")
      ("N_S"     , "number of spins")
      ("N_H"     , "number of holons")
      ("ES"      , "exp(i*pi*s)")
      ;

   C_A = SiteOperator(Basis, QN(0.5,0.5), LatticeCommute::Fermionic);
   C_B = SiteOperator(Basis, QN(0.5,0.5), LatticeCommute::Fermionic);
   P = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   R = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   N_S = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   N_H = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   S = SiteOperator(Basis, QN(0,1), LatticeCommute::Bosonic);
   Q = SiteOperator(Basis, QN(1,0), LatticeCommute::Bosonic);
   I = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   ES = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);

   C_A("spinon", "holon")  = std::sqrt(2.0);
   C_A("holon",  "spinon") = std::sqrt(2.0);

   C_B("spinon", "holon")  =  std::sqrt(2.0);
   C_B("holon",  "spinon") = -std::sqrt(2.0);

   CH_A = adjoint(C_A);
   CH_B = adjoint(C_B);

   P("spinon", "spinon") = -1;
   P("holon",  "holon")  = 1;

   SiteOperator CP_A = prod(C_A, P, QN(0.5,0.5));
   SiteOperator CP_B = prod(C_B, P, QN(0.5,0.5));
   SiteOperator CHP_A = prod(CH_A, P, QN(0.5,0.5));
   SiteOperator CHP_B = prod(CH_B, P, QN(0.5,0.5));

   ES("spinon", "spinon") = std::complex<double>(0, 1);
   ES("holon", "holon") = 1;

   R("spinon", "spinon") = 1;
   R("holon",  "holon")  = 1;

   N_S("spinon", "spinon") = 1;

   N_H("holon", "holon") = 1;

   // this choice of signs makes the 0-component from the Wigner Eckart theorem 
   // equal to the usual definition of Sz.  Note that the "hermitian conjugate" of
   // S is S^\dagger = -S.  This arises from the signs picked up in taking the
   // transpose of a spin 1 operator.
   S("spinon", "spinon") = std::sqrt(0.75);

   Q("holon", "holon")   = std::sqrt(0.75);

   // identity
   I("spinon", "spinon") =  1;
   I("holon",  "holon")  =  1;
   
   Site["I"] = I;
   Site["P"] = P;
   Site["N_S"] = N_S;
   Site["N_H"] = N_H;
   Site[Sym1] = Q;
   Site[Sym2] = S;
   Site["C_A"] = C_A;
   Site["C_B"] = C_B;
   Site["CH_A"] = CH_A;
   Site["CH_B"] = CH_B;
   Site["ES"] = ES;

   Site.set_operator_descriptions(OpDescriptions);

   return Site;
}

inline
LatticeSite FermionSO4_A(std::string const& Sym1 = "Q", std::string const& Sym2 = "S")
{
   LatticeSite Site = FermionSO4_Common(Sym1, Sym2);
   Site["C"] = Site["C_A"];
   Site["C"].set_description("Annihilation operator");
   Site["CH"] = Site["CH_A"];
   Site["CH"].set_description("Creation operator");
   return Site;
}

inline
LatticeSite FermionSO4_B(std::string const& Sym1 = "Q", std::string const& Sym2 = "S")
{
   LatticeSite Site = FermionSO4_Common(Sym1, Sym2);
   Site["C"] = Site["C_B"];
   Site["C"].set_description("Annihilation operator");
   Site["CH"] = Site["CH_B"];
   Site["CH"].set_description("Creation operator");
   return Site;
}

#endif
