// -*- C++ -*- $Id$

#include "lattice/latticesite.h"
#include "quantumnumbers/u1.h"




inline
LatticeSite CreateU1HardcoreBoson(std::string const& Sym = "N")
{
   SymmetryList Symmetry(Sym+":U(1)");
   QuantumNumbers::QNConstructor<QuantumNumbers::U1> QN(Symmetry);
   SiteBasis Basis(Symmetry);
   SiteOperator BH, B, N, P, R, I;
   LatticeSite Site;

   Basis.push_back("empty", QN(0));
   Basis.push_back("boson", QN(1));

   B = SiteOperator(Basis, QN(-1), LatticeCommute::Bosonic);

   B("empty", "boson") = 1.0;
   BH = adjoint(B);

   N = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   P = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   R = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   I = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);

   N = prod(BH, B, QN(0));

   I("empty", "empty") = 1.0;
   I("boson", "boson") = 1.0;

   R = I;

   P("empty", "empty") = 1.0;
   P("boson", "boson") = 1.0;
   
   Site["I"] = I;
   Site["P"] = P;
   Site["R"] = R;
   Site["N"] = N;
   Site["B"] = B;
   Site["BH"] = BH;
   return Site;
}
