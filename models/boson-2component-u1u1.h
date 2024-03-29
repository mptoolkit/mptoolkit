// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/boson-2component-u1u1.h
//
// Copyright (C) 2012-2023 Ian McCulloch <ian@qusim.net>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

#include "lattice/latticesite.h"
#include "quantumnumbers/u1.h"
#include "quantumnumbers/su2.h"




void SetMatElement(SiteOperator& s, int n1, int n2, double x)
{
   std::string q1 = boost::lexical_cast<std::string>(n1);
   std::string q2 = boost::lexical_cast<std::string>(n2);

   int l1 = s.Basis1().LookupOrNeg(q1);
   int l2 = s.Basis2().LookupOrNeg(q2);

   if (l1 >= 0 && l2 >= 0)
      s(q1, q2) = x;
}

std::string Coord(int i, int j)
{
   std::ostringstream S;
   S << '(' << i << ',' << j << ')';
   return S.str();
}

inline
LatticeSite Boson2ComponentU1U1(int MaxN, std::string const& Sym1 = "NA", std::string const& Sym2 = "NB")
{
   SymmetryList Symmetry(Sym1+":U(1),"+Sym2+":U(1)");
   QuantumNumbers::QNConstructor<QuantumNumbers::U1, QuantumNumbers::U1> QN(Symmetry);
   SiteBasis Basis(Symmetry);
   SiteOperator B_A, BH_A, N_A, N2_A;
   SiteOperator B_B, BH_B, N_B, N2_B;
   SiteOperator P, R, Q, I;
   LatticeSite Site("2-component bosons, maximum number of particles per site = " + boost::lexical_cast<std::string>(MaxN));

   OperatorDescriptions OpDescriptions;
   OpDescriptions.add_operators()
      ("I"   , "identity")
      ("R"   , "reflection (trivial)")
      ("B_A" , "annihilation, species A")
      ("BH_A", "creation, species A")
      ("N_A" , "particle number, species A")
      ("N2_A", "Coulomb repulstion N*(N-1), species A")
      ("B_B" , "annihilation, species B")
      ("BH_B", "creation, species B")
      ("N_B" , "particle number, species B")
      ("N2_B", "Coulomb repulstion N*(N-1), species B")

      ;

   // Setup the site basis
   for (int n = 0; n <= MaxN; ++n)
   {
      for (int m = 0; m <= MaxN; ++m)
      {
         std::string q = Coord(n,m);
         Basis.push_back(q, QN(n, m));
      }
   }

   BH_A = SiteOperator(Basis, QN(1, 0), LatticeCommute::Bosonic);
   for (int n = 0; n < MaxN; ++n)
   {
      for (int m = 0; m <= MaxN; ++m)
      {
         std::string q1 = Coord(n+1,m);
         std::string q2 = Coord(n,m);

         int l1 = BH_A.Basis1().LookupOrNeg(q1);
         int l2 = BH_A.Basis2().LookupOrNeg(q2);

         if (l1 >= 0 && l2 >= 0) {
            BH_A(q1, q2) = std::sqrt(double(n + 1));
         } else {
            std::cout << "A: " << n << " " << m << " " << q1 << " " << q2 << '\n';
         }
      }
   }

   BH_B = SiteOperator(Basis, QN(0, 1), LatticeCommute::Bosonic);
   for (int n = 0; n <= MaxN; ++n)
   {
      for (int m = 0; m < MaxN; ++m)
      {
         std::string q1 = Coord(n,m+1);
         std::string q2 = Coord(n,m);

         int l1 = BH_B.Basis1().LookupOrNeg(q1);
         int l2 = BH_B.Basis2().LookupOrNeg(q2);

         if (l1 >= 0 && l2 >= 0) {
            BH_B(q1, q2) = std::sqrt(double(m + 1));
         } else {
            std::cout << "B: " << n << " " << m << " " << q1 << " " << q2 << '\n';
         }
      }
   }

   I = SiteOperator::Identity(Basis);
   Site["I"] = I;
   R = I;
   Site["R"] = R;

   Site["BH_A"] = BH_A;
   Site["BH_B"] = BH_B;

   B_A = adjoint(BH_A);
   Site["B_A"] = B_A;

   B_B = adjoint(BH_B);
   Site["B_B"] = B_B;

   N_A = prod(BH_A, B_A, QN(0,0));
   Site["N_A"] = N_A;
   N2_A = prod(N_A, N_A-I, QN(0,0));
   Site["N2_A"] = N2_A;

   N_B = prod(BH_B, B_B, QN(0,0));
   Site["N_B"] = N_B;
   N2_B = prod(N_B, N_B-I, QN(0,0));
   Site["N2_B"] = N2_B;

   Site.set_operator_descriptions(OpDescriptions);

   return Site;
}
