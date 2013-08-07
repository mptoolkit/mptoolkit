// -*- C++ -*- $Id$

// Bose-Hubbard model with two species of bosons,
// symmetric under interchange of species.

#include "siteoperator/siteoperator.h"
#include "quantumnumbers/u1.h"
#include "quantumnumbers/z2.h"
#include "siteoperator/block.h"

typedef Block<SiteOperator> SiteBlock;

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
SiteBlock CreateBoseHubbard2BosonsU1Z2Site(int MaxN, std::string const& Sym1 = "N", std::string const& Sym2 = "Z")
{
   SymmetryList Symmetry(Sym1+":U(1),"+Sym2+":Z2");
   QuantumNumbers::QNConstructor<QuantumNumbers::U1, QuantumNumbers::Z2> QN(Symmetry);
   SiteBasis Basis(Symmetry);
   SiteOperator BH_A, BH_S, I, Z;
   SiteBlock Block;

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
   Block["I"] = I;
   Block["R"] = I;
   Block["P"] = I;

   Block["Z"]     = Z;
   Block["BH_A"]  = BH_A;
   Block["BH_S"]  = BH_S;
   Block["B_A"]   = adjoint(BH_A);
   Block["B_S"]   = adjoint(BH_S);

   Block["BH2_A"] = prod(Block["BH_A"], Block["BH_A"], QN(2, 1));;
   Block["B2_A"]  = prod(Block["B_A"], Block["B_A"], QN(-2, 1));
   Block["BH2_S"] = prod(Block["BH_S"], Block["BH_S"], QN(2, 1));;
   Block["B2_S"]  = prod(Block["B_S"], Block["B_S"], QN(-2, 1));
   Block["N_A"]   = prod(Block["BH_A"], Block["B_A"], QN(0,1));
   Block["N2_A"]  = prod(Block["N_A"], Block["N_A"]-Block["I"], QN(0,1));
   Block["N_S"]   = prod(Block["BH_S"], Block["B_S"], QN(0,1));
   Block["N2_S"]  = prod(Block["N_S"], Block["N_S"]-Block["I"], QN(0,1));
   Block["N"]     = Block["N_A"] + Block["N_S"];
   
   DEBUG_TRACE(BH_A)(B_A)(BH_S)(B_S)(I)(N_A)(N2_A)(N_S)(N2_S);

   return Block;
}
