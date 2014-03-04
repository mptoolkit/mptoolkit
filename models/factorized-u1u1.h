// -*- C++ -*- $Id$
//
// factorized-u1u1.h
//
// Site blocks for bosonic and fermionic models factorized into 
// (0 up) \otimes (0 down) sites.
//

#include "lattice/latticesite.h"
#include "quantumnumbers/u1.h"




inline
LatticeSite FactorizedSite(int MyN, half_int MySz, 
                         std::string const& Sym1 = "N", 
                         std::string const& Sym2 = "Sz")
{
   SymmetryList Symmetry(Sym1+":U(1),"+Sym2+":U(1)");
   QuantumNumbers::QNConstructor<QuantumNumbers::U1,QuantumNumbers::U1> QN(Symmetry);
   SiteBasis Basis(Symmetry);
   SiteOperator CH, C, P, R, N, I, Sz;
   LatticeSite Site;

   Basis.push_back("0", QN(0,0));
   Basis.push_back("1", QN(MyN,MySz));

   CH = SiteOperator(Basis, QN(MyN,MySz), LatticeCommute(LatticeCommute::Fermionic) * MyN);
   C = SiteOperator(Basis, QN(-MyN,-MySz), LatticeCommute(LatticeCommute::Fermionic) * MyN);
   P = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   R = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   N = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   Sz = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   I = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);

   // annihilate particle
   C("0", "1") = std::sqrt(double(MyN));
   
   // create particle
   CH = adjoint(C);

   // identity
   I("0", "0") = 1;
   I("1", "1") = 1;

   // parity = (-1)^N
   P("0", "0") = 1;
   P("1", "1") = minus1pow(MyN);

   R = I;

   // particle number  
   N = prod(CH, C, QN(0,0));

   // z-component of spin
   Sz("1", "1") = MySz.to_double();

   Site["I"] = I;
   Site["P"] = P;
   Site["R"] = R;
   Site["N"] = N;
   Site["Sz"] = Sz;
   Site["CH"] = CH;
   Site["C"] = C;
   return Site;
}
