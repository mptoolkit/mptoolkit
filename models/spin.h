// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/spin.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

inline
LatticeSite SpinSite(half_int Spin)
{
   SymmetryList Symmetry("S:Null");
   SiteBasis Basis(Symmetry);
   SiteOperator Sp, Sm, Sz, Sx, Sy, R, P, I;
   LatticeSite Site("Spin "+to_string_fraction(Spin));

   QuantumNumbers::QuantumNumber q(Symmetry); // no symmetries, only one quantum number

   std::map<half_int, std::string> SpinBasis;
   for (half_int s = -Spin; s <= Spin; ++s)
   {
      SpinBasis[s] = boost::lexical_cast<std::string>(s);
      Basis.push_back(SpinBasis[s], q);
   }

   OperatorDescriptions OpDescriptions;
   OpDescriptions.add_operators()
      ("I"   , "identity")
      ("R"   , "reflection")
      ("P"   , "fermion parity")
      ("Sp"  , "raising operator")
      ("Sm"  , "lowering operator")
      ("Sx"  , "x-component of spin")
      ("Sy"  , "y-component of spin")
      ("Sz"  , "z-component of spin")
      ;

   if (Spin == 0.5)
   {
      OpDescriptions.add_operators()
         ("X", "Pauli Sigma-X")
         ("Y", "Pauli Sigma-Y")
         ("Z", "Pauli Sigma-Z")
         ;
   }


   Sp = SiteOperator(Basis, q, LatticeCommute::Bosonic);
   Sm = SiteOperator(Basis, q, LatticeCommute::Bosonic);
   Sx = SiteOperator(Basis, q, LatticeCommute::Bosonic);
   Sy = SiteOperator(Basis, q, LatticeCommute::Bosonic);
   Sz = SiteOperator(Basis, q, LatticeCommute::Bosonic);

   P = SiteOperator(Basis, q, LatticeCommute::Bosonic);
   R = SiteOperator(Basis, q, LatticeCommute::Bosonic);
   I = SiteOperator(Basis, q, LatticeCommute::Bosonic);

   for (half_int s = -Spin; s <= Spin; ++s)
   {
      I(SpinBasis[s], SpinBasis[s]) = 1.0;
      P(SpinBasis[s], SpinBasis[s]) = 1.0;
      R(SpinBasis[s], SpinBasis[s]) = 1.0;
      Sz(SpinBasis[s], SpinBasis[s]) = s.to_double();
   }

   for (half_int s = -Spin; s < Spin; ++s)
   {
      Sp(SpinBasis[s+1], SpinBasis[s]) = std::sqrt((Spin - s) * (Spin + s + 1));
   }

   Sm = adjoint(Sp);

   Sx = 0.5 * (Sp + Sm);
   Sy = std::complex<double>(0.0, -0.5) * (Sp - Sm);

   Site.arg("Spin") = Spin.to_double();

   Site.func("Ux")(arg("theta") = math_const::pi) = "exp(theta*i*Sx)";
   Site.func("Uy")(arg("theta") = math_const::pi) = "exp(theta*i*Sy)";
   Site.func("Uz")(arg("theta") = math_const::pi) = "exp(theta*i*Sz)";

   Site["I"] = I;
   Site["P"] = P;
   Site["R"] = R;
   Site["Sp"] = Sp;
   Site["Sm"] = Sm;
   Site["Sx"] = Sx;
   Site["Sy"] = Sy;
   Site["Sz"] = Sz;
   if (Spin == 0.5)
   {
      Site["X"] = 2*Sx;
      Site["Y"] = 2*Sy;
      Site["Z"] = 2*Sz;
   }

   Site.set_operator_descriptions(OpDescriptions);

   return Site;
}
