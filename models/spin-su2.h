// -*- C++ -*- $Id$

#include "siteoperator/siteoperator.h"
#include "siteoperator/block.h"
#include "quantumnumbers/su2.h"

typedef Block<SiteOperator> SiteBlock;

inline
SiteBlock CreateSU2SpinSite(half_int Spin, std::string const& Sym = "S")
{
   SymmetryList Symmetry(Sym+":SU(2)");
   QuantumNumbers::QNConstructor<QuantumNumbers::SU2> QN(Symmetry);
   SiteBasis Basis(Symmetry);
   SiteOperator S, Sn, R, P, I, Q, T, F;
   SiteBlock Block;

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
   
   Block["I"] = I;
   Block["P"] = P;
   Block["R"] = R;
   Block["S"] = S;
   Block["Sn"] = Sn;
   Block["Q"] = Q;
   Block["T"] = T;
   Block["F"] = F;
   return Block;
}
