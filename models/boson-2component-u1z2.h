// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/boson-2component-u1z2.h
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

// Bose-Hubbard model with two species of bosons,
// symmetric under interchange of species.

#include "lattice/latticesite.h"
#include "quantumnumbers/u1.h"
#include "quantumnumbers/z2.h"

// naming of the basis states, for (i,j) being the number
// of symmetric and antisymmetric bosons respectively
std::string Coord(int i, int j)
{
   std::ostringstream S;
   S << '(' << i << ',' << j << ')';
   return S.str();
}

int Parity(int m)
{
   // FIXME: this is the negative of what it should be
   return (m%2)*2-1;
}

inline
LatticeSite Boson2ComponentU1Z2(int MaxN, std::string const& Sym1 = "N", std::string const& Sym2 = "Z")
{
   SymmetryList Symmetry(Sym1+":U(1),"+Sym2+":Z2");
   QuantumNumbers::QNConstructor<QuantumNumbers::U1, QuantumNumbers::Z2> QN(Symmetry);
   SiteBasis Basis(Symmetry);
   LatticeSite Site("Two-component bosons in the U(1) x Z2 basis");

   SiteOperator& BH_A = Site["BH_A"];
   SiteOperator& BH_S = Site["BH_S"];
   SiteOperator& B_A = Site["B_A"];
   SiteOperator& B_S = Site["B_S"];
   SiteOperator& N_A = Site["N_A"];
   SiteOperator& N_S = Site["N_S"];
   SiteOperator& N = Site["N"];
   SiteOperator& N2_A = Site["N2_A"];
   SiteOperator& N2_S = Site["N2_S"];
   SiteOperator& I = Site["I"];
   SiteOperator& Z = Site["Z"];
   SiteOperator& D = Site["D"];
   SiteOperator& D2 = Site["D2"];  // square of the order parameter

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
   for (int n = 0; n <= MaxN; ++n)
   {
      for (int m = 0; m <= MaxN; ++m)
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
   Site["R"] = I;
   Site["P"] = I;

   B_A = adjoint(BH_A);
   B_S = adjoint(BH_S);

   N_A = BH_A*B_A;
   N_S = BH_S*B_S;

   N = N_A + N_S;
   
   N2_A = N_A * (N_A-I);
   N2_S = N_S * (N_S-I);

   D = BH_S*B_A + BH_A*B_S;
   D2 = D*D;

   return Site;
}
