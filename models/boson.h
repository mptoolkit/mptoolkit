// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/boson.h
//
// Copyright (C) 2014-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#if !defined(MPTOOLKIT_MODELS_BOSON_H)
#define MPTOOLKIT_MODELS_BOSON_H

#include "lattice/latticesite.h"
#include "models/boson-u1.h"    // for DefaultMaxN

void SetMatElement(SiteOperator& s, int n1, int n2, double x)
{
   std::string q1 = boost::lexical_cast<std::string>(n1);
   std::string q2 = boost::lexical_cast<std::string>(n2);

   int l1 = s.Basis1().LookupOrNeg(q1);
   int l2 = s.Basis2().LookupOrNeg(q2);

   if (l1 >= 0 && l2 >= 0)
      s(q1, q2) = x;
}

inline
LatticeSite Boson(int MaxN)
{
   SymmetryList Symmetry("N:Null");
   QuantumNumbers::QuantumNumber QNum(Symmetry); // no symmetries, only one quantum number
   SiteBasis Basis(Symmetry);
   SiteOperator B, BH, P, R, N, N2, Q, I, U;
   LatticeSite Site;

   // Setup the site basis
   for (int n = 0; n <= MaxN; ++n)
   {
      std::string q = boost::lexical_cast<std::string>(n);
      Basis.push_back(q, QNum);
   }

   BH = SiteOperator(Basis, QNum, LatticeCommute::Bosonic);
   for (int n = 0; n <= MaxN; ++n)
   {
     SetMatElement(BH, n+1, n, std::sqrt(double(n + 1)));
   }

   Site["BH"] = BH;

   B = adjoint(BH);
   Site["B"] = B;
   I = SiteOperator::Identity(Basis);
   Site["I"] = I;
   N = prod(BH, B, QNum);
   Site["N"] = N;
   N2 = prod(N, N-I, QNum);
   Site["N2"] = N2;
   R = I;
   Site["R"] = R;

   U = SiteOperator(Basis, QNum, LatticeCommute::Bosonic);
   for (int n = 0; n <= MaxN; ++n)
   {
      SetMatElement(U, n, n, minus1pow(n));
   }
   Site["U"] = U;

   // projections onto each state
   for (int n = 0; n <= MaxN; ++n)
   {
      SiteOperator X(Basis, QNum, LatticeCommute::Bosonic);
      SetMatElement(X, n, n, 1);
      std::string OpName = std::string("P_")
        + boost::lexical_cast<std::string>(n);
      Site[OpName] = X;
   }

   DEBUG_TRACE(BH)(B)(I)(N)(N2)(U);

   return Site;
}

#endif
