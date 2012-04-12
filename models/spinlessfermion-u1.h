// -*- C++ -*- $Id$

#include "siteoperator/siteoperator.h"
#include "quantumnumbers/u1.h"
#include "siteoperator/block.h"

typedef Block<SiteOperator> SiteBlock;

inline
SiteBlock CreateU1SpinlessFermion(std::string const& Sym = "N")
{
   SymmetryList Symmetry(Sym+":U(1)");
   QuantumNumbers::QNConstructor<QuantumNumbers::U1> QN(Symmetry);
   SiteBasis Basis(Symmetry);
   SiteOperator CH, C, N, P, R, I;
   SiteBlock Block;

   Basis.push_back("empty", QN(0));
   Basis.push_back("fermion", QN(1));

   C = SiteOperator(Basis, QN(-1), LatticeCommute::Fermionic);

   C("empty", "fermion") = 1.0;
   CH = adjoint(C);

   N = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   P = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   R = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   I = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);

   N = prod(CH, C, QN(0));

   I("empty", "empty") = 1.0;
   I("fermion", "fermion") = 1.0;

   R = I;

   P("empty", "empty") = 1.0;
   P("fermion", "fermion") = -1.0;
   
   Block["I"] = I;
   Block["P"] = P;
   Block["R"] = R;
   Block["N"] = N;
   Block["C"] = C;
   Block["CH"] = CH;
   Block["CP"] = prod(C, P, QN(-1));
   Block["CHP"] = prod(CH, P, QN(1));
   return Block;
}
