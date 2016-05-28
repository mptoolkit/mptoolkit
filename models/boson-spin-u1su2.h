// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/boson-spin-u1su2.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Reseach publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

#include "lattice/latticesite.h"
#include "quantumnumbers/u1.h"
#include "quantumnumbers/su2.h"




void SetMatElement(SiteOperator& s, int n1, int s1, int n2, int s2, double x)
{
   std::string q1 = boost::lexical_cast<std::string>(n1) 
      + "," + boost::lexical_cast<std::string>(s1);
   std::string q2 = boost::lexical_cast<std::string>(n2) 
      + "," + boost::lexical_cast<std::string>(s2);
   
   int l1 = s.Basis1().LookupOrNeg(q1);
   int l2 = s.Basis2().LookupOrNeg(q2);

   if (l1 >= 0 && l2 >= 0)
      s(q1, q2) = x;
}

inline
LatticeSite CreateU1SU2BoseHubbardSite(int MaxN, int MaxS, 
				     std::string const& Sym1 = "N", 
				     std::string const& Sym2 = "S")
{
   // The matrix elements for BH are only defined up to N=5
   CHECK(MaxN <= 5 && MaxS <= 5); 
   
   SymmetryList Symmetry(Sym1+":U(1),"+Sym2+":SU(2)");
   QuantumNumbers::QNConstructor<QuantumNumbers::U1,QuantumNumbers::SU2> QN(Symmetry);
   SiteBasis Basis(Symmetry);
   SiteOperator B, BH, P, R, S, S2, N, N2, Q, I;
   LatticeSite Site;

   // Setup the site basis
   for (int n = 0; n <= MaxN; ++n)
   {
      for (int s = n%2; s <= std::min(n, MaxS); s += 2)
      {
	 std::string q = boost::lexical_cast<std::string>(n) 
	    + "," + boost::lexical_cast<std::string>(s);
         Basis.push_back(q, QN(n,s));
      }
   }

   BH = SiteOperator(Basis, QN(1,1), LatticeCommute::Bosonic);
   for (int n = 0; n <= MaxN; ++n)
   {
      for (int s = n%2; s <= std::min(n, MaxS); s += 2)
      {
         SetMatElement(BH, n+1, s+1, n, s, std::sqrt(double(n + 1)));
         SetMatElement(BH, n+1, s-1, n, s, std::sqrt(double(n + 1)));
      }
   }

   // kludge the spin factors
   SetMatElement(BH, 3,1, 2,0, std::sqrt(3.0 * 5.0 / 9.0));
   SetMatElement(BH, 3,1, 2,2, std::sqrt(3.0 * 4.0 / 9.0));

   SetMatElement(BH, 4,2, 3,1, std::sqrt(4.0 * 7.0 / 10.0));
   SetMatElement(BH, 4,2, 3,3, std::sqrt(4.0 * 3.0 / 10.0));

   SetMatElement(BH, 5,1, 4,0, std::sqrt(5.0 * 7.0 / 15.0));
   SetMatElement(BH, 5,1, 4,2, std::sqrt(5.0 * 8.0 / 15.0));

   SetMatElement(BH, 5,3, 4,2, std::sqrt(5.0 * 27.0 / 35.0));
   SetMatElement(BH, 5,3, 4,4, std::sqrt(5.0 * 8.0 / 35.0));

   Site["BH"] = BH;

   B = adjoint(BH);
   Site["B"] = B;
   I = SiteOperator::Identity(Basis);
   Site["I"] = I;
   N = sqrt(3.0) * prod(BH, B, QN(0,0));
   Site["N"] = N;
   N2 = prod(N, N-I, QN(0,0));
   Site["N2"] = N2;
   R = I;
   Site["R"] = R;
   S = sqrt(2.0) * prod(BH, B, QN(0,1));
   Site["S"] = S;
   S2 = -sqrt(3.0) * prod(S, S, QN(0,0));
   Site["S2"] = S2;
   Q = prod(S, S, QN(0,2));
   Site["Q"] = Q;

   // projections onto each state
   for (int n = 0; n <= MaxN; ++n)
   {
      for (int s = n%2; s <= std::min(n, MaxS); s += 2)
      {
         SiteOperator X(Basis, QN(0,0), LatticeCommute::Bosonic);
         SetMatElement(X, n, s, n, s, 1);
	 std::string OpName = std::string("P(")
	    + boost::lexical_cast<std::string>(n) + "," 
	    + boost::lexical_cast<std::string>(s) + ")";
	 Site[OpName] = X;
      }
   }

   DEBUG_TRACE(BH)(B)(I)(N)(N2)(S)(S2)(Q);

   return Site;
}
