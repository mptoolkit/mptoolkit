// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/contrib/o2-u1.h
//
// Copyright (C) 2004-2022 Ian McCulloch <ian@qusim.net>
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
LatticeSite O2U1(int Spin, std::string const& Sym = "Sz")
{
   SymmetryList Symmetry(Sym+":U(1)");
   QuantumNumbers::QNConstructor<QuantumNumbers::U1> QN(Symmetry);
   SiteBasis Basis(Symmetry);
   SiteOperator Up, Um, Sz, R, P, I;
   LatticeSite Site("U(1) O(2) "+to_string_fraction(Spin));

   std::map<int, std::string> SpinBasis;
   for (int s = -Spin; s <= Spin; ++s)
   {
      SpinBasis[s] = boost::lexical_cast<std::string>(s);
      Basis.push_back(SpinBasis[s], QN(s));
   }

   OperatorDescriptions OpDescriptions;
   OpDescriptions.add_operators()
      ("I"   , "identity")
      ("R"   , "reflection")
      ("P"   , "fermion parity")
      ("Up"  , "raising operator")
      ("Um"  , "lowering operator")
      ("Sz"  , "z-component of spin")
      ;

   Up = SiteOperator(Basis, QN(1), LatticeCommute::Bosonic);
   Um = SiteOperator(Basis, QN(-1), LatticeCommute::Bosonic);
   Sz = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);

   P = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   R = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);
   I = SiteOperator(Basis, QN(0), LatticeCommute::Bosonic);

   for (int s = -Spin; s <= Spin; ++s)
   {
      I(SpinBasis[s], SpinBasis[s]) = 1.0;
      P(SpinBasis[s], SpinBasis[s]) = 1.0;
      R(SpinBasis[s], SpinBasis[s]) = 1.0;
      Sz(SpinBasis[s], SpinBasis[s]) = s;
   }

   // Sp and Sm operators
   for (int s = -Spin; s < Spin; ++s)
   {
      Up(SpinBasis[s+1], SpinBasis[s]) = 1.0;
   }
   Um = adjoint(Up);

   // Example of defining a named constant (argument)
   Site.arg("Cutoff") = Spin;

   // Example of defining a function.  The first parameter has a default value
   Site.func("Uz")(arg("theta") = math_const::pi) = "exp(theta*i*Sz)";

   Site["I"] = I;
   Site["P"] = P;
   Site["R"] = R;
   Site["Up"] = Up;
   Site["Um"] = Um;
   Site["Sz"] = Sz;

   Site.set_operator_descriptions(OpDescriptions);

   return Site;
}
