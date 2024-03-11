// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/contrib/spin-su2u1.h
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
//
#include "lattice/latticesite.h"
#include "quantumnumbers/su2.h"
#include "quantumnumbers/u1.h"

// SU(3) spin chain in the fundamental representation,
// broken down to SU(1)xU(1)
//
// Start from the Gell-Mann matrices
//
// L1 = [  0  1  0 ]
//      [  1  0  0 ]
//      [  0  0  0 ]
//
// L2 = [  0 -i  0 ]
//      [  i  0  0 ]
//      [  0  0  0 ]
//
// L3 = [  1  0  0 ]
//      [  0 -1  0 ]
//      [  0  0  0 ]
//
// L4 = [  0  0  1 ]
//      [  0  0  0 ]
//      [  1  0  0 ]
//
// L5 = [  0  0 -i ]
//      [  0  0  0 ]
//      [  i  0  0 ]
//
// L6 = [  0  0  0 ]
//      [  0  0  1 ]
//      [  0  1  0 ]
//
// L7 = [  0  0  0 ]
//      [  0  0 -i ]
//      [  0  i  0 ]
//
// L8 = (1/sqrt(3)) [  1  0  0 ]
//                  [  0  1  0 ]
//                  [  0  0 -2 ]
//
// We take (L1,L2,L3) to be SU(2) generators, via
//
// Sx = L1/2
// Sy = L2/2
// Sz = L3/2
//
// and sqrt(3) * L8 as the U(1) generator.  Note that L8 commutes with Sx,Sy,Sz, so these symmetries are independent.
// The basis is an SU(2) doublet and an SU(2) singlet, that have U(1) charges 1 and -2 respectively.
//
// The Hamiltonian
// ===============
//
// The 8 generators split into
//
// inner(S,S) = L1.L1 + L2.L2 + L3.L3
// inner(U,U) + inner(V,V) = L4.L4 + L5.L5 + L6.L6 + L7.L7
// L8.L8

inline
LatticeSite SpinSU2U1()
{
   SymmetryList Symmetry("S:SU(2),Qz:U(1)");
   QuantumNumbers::QNConstructor<QuantumNumbers::SU2, QuantumNumbers::U1>
      QN(Symmetry);
   SiteBasis Basis(Symmetry);
   SiteOperator L8;
   SiteOperator S, U, V, Qz;

   SiteOperator R, P, I;
   LatticeSite Site;

   // The quantum numbers of the basis are L3/2 and sqrt(3)*L8
   Basis.push_back("s", QN(0.5, 1));
   Basis.push_back("0", QN(0, -2));

   OperatorDescriptions OpDescriptions;
   OpDescriptions.add_operators()
      ("I"   , "identity")
      ("R"   , "reflection")
      ("P"   , "fermion parity")
      ("L8"  , "U(1) generator lambda_8")
      ("Qz"  , "U(1) generator scaled to have integer eigenvalues 1,1,-2")
      ("S"   , "SU(2) spin operator")
      ("U"   , "Spin creation operator")
      ("V"   , "Spin annihilation operator")
      ;

   L8 = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);

   S = SiteOperator(Basis, QN(1,0), LatticeCommute::Bosonic);

   Qz = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);

   // Annihilate a spin to a (0,-2) state
   V = SiteOperator(Basis, QN(0.5,-3), LatticeCommute::Bosonic);

   // Create a spin-1/2, conjugate of V
   U = SiteOperator(Basis, QN(0.5,3), LatticeCommute::Bosonic);

   P = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   R = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   I = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);

   S("s", "s")    = std::sqrt(3.0);

   L8("s", "s")   =  1.0 / std::sqrt(3.0);
   L8("0", "0")   = -2.0 / std::sqrt(3.0);

   Qz("s", "s")   =  1.0;
   Qz("0", "0")   = -2.0;

   V("0", "s")    =  2.0;
   U = adjoint(V);

   I("s", "s")    =  1.0;
   I("0", "0")    =  1.0;

   P = I;
   R = I;

   Site["I"] = I;
   Site["P"] = P;
   Site["R"] = R;
   Site["L8"] = L8;
   Site["Qz"] = Qz;
   Site["S"] = S;
   Site["U"] = U;
   Site["V"] = V;

   Site.set_operator_descriptions(OpDescriptions);

   return Site;
}
