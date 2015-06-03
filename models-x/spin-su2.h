// -*- C++ -*- $Id$

#include "lattice/latticesite.h"
#include "quantumnumbers/su2.h"

inline
LatticeSite
CreateSU2SpinSite(half_int Spin, std::string const& Sym = "S")
{
   SymmetryList Symmetry(Sym+":SU(2)");
   QuantumNumbers::QNConstructor<QuantumNumbers::SU2> QN(Symmetry);
   SiteBasis Basis(Symmetry);
   SiteOperator S, Sn, R, P, I, Q, T, F;
   LatticeSite Site;

   std::string SpinName = boost::lexical_cast<std::string>(Spin);
   Basis.push_back(SpinName, QN(Spin));

   S = SiteOperator(Basis, QN(1), LatticeCommute::Bosonic);
   P = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   R = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   I = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   Sn = SiteOperator(Basis, QN(1), LatticeCommute::Bosonic);

   S(SpinName, SpinName) = sqrt(Spin * (Spin+1));
   Sn(SpinName, SpinName) = 1;
   I(SpinName, SpinName) = 1.0;
   P = I;
   R = I;

   Q = prod(S, S, QN(2));
   T = prod(Q, S, QN(3));
   F = prod(Q, Q, QN(4));
   
   Site["I"] = I;
   Site["P"] = P;
   Site["R"] = R;
   Site["S"] = S;
   Site["Sn"] = Sn;
   Site["Q"] = Q;
   Site["T"] = T;
   Site["F"] = F;
   return Site;
}
