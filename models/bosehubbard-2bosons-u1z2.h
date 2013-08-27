// -*- C++ -*- $Id$

// Bose-Hubbard model with two species of bosons,
// symmetric under interchange of species.

#include "siteoperator/latticesite.h"
#include "quantumnumbers/u1.h"
#include "quantumnumbers/z2.h"




std::string Coord(int i, int j)
{
   std::ostringstream S;
   S << '(' << i << ',' << j << ')';
   return S.str();
}

int Parity(int m)
{
   return (m%2)*2-1;
}

inline
LatticeSite CreateBoseHubbard2BosonsU1Z2Site(int MaxN, std::string const& Sym1 = "N", std::string const& Sym2 = "Z")
{
   SymmetryList Symmetry(Sym1+":U(1),"+Sym2+":Z2");
   QuantumNumbers::QNConstructor<QuantumNumbers::U1, QuantumNumbers::Z2> QN(Symmetry);
   SiteBasis Basis(Symmetry);
   SiteOperator BH_A, BH_S, I, Z;
   LatticeSite Site;

   // Setup the site basis
   for (int n = 0; n <= MaxN; ++n)
   {
      for (int m = 0; m <= MaxN; ++m)
      {
	 Basis.push_back(Coord(n, m), QN(n+m, Parity(m)));
      }
   }

   // creation operators in the anti-symmetric and symmetric channels
   BH_A = SiteOperator(Basis, QN(1,-1), LatticeCommute::Bosonic);
   BH_S = SiteOperator(Basis, QN(1, 1), LatticeCommute::Bosonic);
   // parity
   Z = SiteOperator(Basis, QN(0, 1), LatticeCommute::Bosonic);
   for (int n = 0; n < MaxN; ++n)
   {
      for (int m = 0; m < MaxN; ++m)
      {
	 int l1 = Basis.LookupOrNeg(Coord(n, m+1));
	 int l2 = Basis.LookupOrNeg(Coord(n, m));
	 if (l1 >= 0 && l2 >= 0) 
	 {
	    BH_A(Coord(n, m+1), Coord(n, m)) = std::sqrt(double(m + 1));
	 }
	 
	 l1 = Basis.LookupOrNeg(Coord(n+1, m));
	 l2 = Basis.LookupOrNeg(Coord(n, m));
	 if (l1 >= 0 && l2 >= 0) 
	 {
	    BH_S(Coord(n+1, m), Coord(n, m)) = std::sqrt(double(n + 1));
	 }

	 Z(Basis.LookupOrNeg(Coord(n, m)), Basis.LookupOrNeg(Coord(n, m))) = Parity(m);
      }
   }

   I = SiteOperator::Identity(Basis);
   Site["I"] = I;
   Site["R"] = I;
   Site["P"] = I;

   Site["Z"]     = Z;
   Site["BH_A"]  = BH_A;
   Site["BH_S"]  = BH_S;
   Site["B_A"]   = adjoint(BH_A);
   Site["B_S"]   = adjoint(BH_S);

   Site["BH2_A"] = prod(Site["BH_A"], Site["BH_A"], QN(2, 1));;
   Site["B2_A"]  = prod(Site["B_A"], Site["B_A"], QN(-2, 1));
   Site["BH2_S"] = prod(Site["BH_S"], Site["BH_S"], QN(2, 1));;
   Site["B2_S"]  = prod(Site["B_S"], Site["B_S"], QN(-2, 1));
   Site["N_A"]   = prod(Site["BH_A"], Site["B_A"], QN(0,1));
   Site["N2_A"]  = prod(Site["N_A"], Site["N_A"]-Site["I"], QN(0,1));
   Site["N_S"]   = prod(Site["BH_S"], Site["B_S"], QN(0,1));
   Site["N2_S"]  = prod(Site["N_S"], Site["N_S"]-Site["I"], QN(0,1));
   Site["N"]     = Site["N_A"] + Site["N_S"];
   
   DEBUG_TRACE(Z)(BH_A)(BH_S);

   return Site;
}
