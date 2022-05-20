// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/spin-u1u1.h
//
// Copyright (C) 2004-2022 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
#include "quantumnumbers/u1.h"

// SU(4) spin chain in the fundamental representation,
// broken down to U(1)xU(1)xU(1)
//
// This uses the basis from http://www.ejtp.com/articles/ejtpv10i28.pdf
//
// L1 = [  0  1  0  0 ]
//      [  1  0  0  0 ]
//      [  0  0  0  0 ]
//      [  0  0  0  0 ]
//
// L2 = [  0 -i  0  0 ]
//      [  i  0  0  0 ]
//      [  0  0  0  0 ]
//      [  0  0  0  0 ]
//
// L3 = [  1  0  0  0 ]
//      [  0 -1  0  0 ]
//      [  0  0  0  0 ]
//      [  0  0  0  0 ]
//
// L4 = [  0  0  1  0 ]
//      [  0  0  0  0 ]
//      [  1  0  0  0 ]
//      [  0  0  0  0 ]
//
// L5 = [  0  0 -i  0 ]
//      [  0  0  0  0 ]
//      [  i  0  0  0 ]
//      [  0  0  0  0 ]
//
// L6 = [  0  0  0  0 ]
//      [  0  0  1  0 ]
//      [  0  1  0  0 ]
//      [  0  0  0  0 ]
//
// L7 = [  0  0  0  0 ]
//      [  0  0 -i  0 ]
//      [  0  i  0  0 ]
//      [  0  0  0  0 ]
//
// L8 = (1/sqrt(3)) [  1  0  0  0 ]
//                  [  0  1  0  0 ]
//                  [  0  0 -2  0 ]
//                  [  0  0  0  0 ]
//
// L9 = [  0  0  0  1 ]
//      [  0  0  0  0 ]
//      [  0  0  0  0 ]
//      [  1  0  0  0 ]
//
// L10 = [  0  0  0 -i ]
//       [  0  0  0  0 ]
//       [  0  0  0  0 ]
//       [  i  0  0  0 ]
//
// L11 = [  0  0  0  0 ]
//       [  0  0  0  1 ]
//       [  0  0  0  0 ]
//       [  0  1  0  0 ]
//
// L12 = [  0  0  0  0 ]
//       [  0  0  0 -i ]
//       [  0  0  0  0 ]
//       [  0  i  0  0 ]
//
// L13 = [  0  0  0  0 ]
//       [  0  0  0  0 ]
//       [  0  0  0  1 ]
//       [  0  0  1  0 ]
//
// L14 = [  0  0  0  0 ]
//       [  0  0  0  0 ]
//       [  0  0  0 -i ]
//       [  0  0  i  0 ]
//
// L15 = (1/sqrt(6)) [  1  0  0  0 ]
//                   [  0  1  0  0 ]
//                   [  0  0  1  0 ]
//                   [  0  0  0 -3 ]
//
// We take the U(1) labels to be L3, sqrt(3)*L8, and sqrt(6)*L15.
//
// A conventient representation is
// Tp = L1 + i*L2
// Tm = L1 - i*L2
// Vp = L4 + i*L5
// Vp = L4 - i*L5
// Up = L6 + i*L7
// Um = L6 - i*L7
// Wp = L9 + i*L10
// Wm = L9 - i*L10
// Xp = L11 + i*L12
// Xm = L11 - i*L12
// Yp = L13 + i*L14
// Ym = L13 - i*L14
//
// We can assemble these into SU(2) generators by adding:
// Tz = L3
// Vz =
//
// These split into SU(2)xSU(2) in 3 different ways:
// [T,Y] = 0
// [V,X] = 0
// [U,W] = 0
//

inline
LatticeSite SpinU1U1U1()
{
   SymmetryList Symmetry("L3:U(1),L8:U(1),L15:U(1)");
   QuantumNumbers::QNConstructor<QuantumNumbers::U1, QuantumNumbers::U1, QuantumNumbers::U1>
      QN(Symmetry);
   SiteBasis Basis(Symmetry);
   SiteOperator L3, L8, L15;
   SiteOperator Tp, Tm, Up, Um, Vp, Vm, Wp, Wm, Xp, Xm, yp, Ym;

   SiteOperator R, P, I, U;
   LatticeSite Site;

   Basis.push_back("0", QN(1, 1, 1));
   Basis.push_back("1", QN(-1, 1, 1));
   Basis.push_back("2", QN(0, -2, 1));
   Basis.push_back("3", QN(0, 0, -3));

   OperatorDescriptions OpDescriptions;
   OpDescriptions.add_operators()
      ("I"   , "identity")
      ("R"   , "reflection")
      ("P"   , "fermion parity")
      ("U"   , "unitary rotation in lambda_8 direction")
      ("L3"  , "U(1) generator lambda_3")
      ("L8"  , "U(1) generator lambda_8")
      ("L15" , "U(1) generator lambda_15")
      ("M8"  , "U(1) generator lambda_3 scaled to have integer eigenvalues 1,-1,0,0")
      ("M8"  , "U(1) generator lambda_8 scaled to have integer eigenvalues 1,1,-2,0")
      ("M15" , "U(1) generator lambda_15 scaled to have integer eigenvalues 1,1,1,-3")
      ("Tp"  , "T raising operator")
      ("Tm"  , "T lowering operator")
      ("Vp"  , "V raising operator")
      ("Vm"  , "V lowering operator")
      ("Up"  , "U raising operator")
      ("Um"  , "U lowering operator")
      ("Wp"  , "W raising operator")
      ("Wm"  , "W lowering operator")
      ("Xp"  , "X raising operator")
      ("Xm"  , "X lowering operator")
      ("Yp"  , "Y raising operator")
      ("Ym"  , "Y lowering operator")
      ;

   L3 = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   L8 = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);

   Sz = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   Qz = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);

   // T+ = L1 + i*L2
   // T- = L1 - i*L2
   Tp = SiteOperator(Basis, QN(1,0), LatticeCommute::Bosonic);
   Tm = SiteOperator(Basis, QN(-1,0), LatticeCommute::Bosonic);

   // V+ = L4 + i*L5
   // V- = L4 - i*L5
   Vp = SiteOperator(Basis, QN(0.5,3), LatticeCommute::Bosonic);
   Vm = SiteOperator(Basis, QN(-0.5,-3), LatticeCommute::Bosonic);

   // U+ = L6 + i*L7
   // U- = L6 - i*L7
   Up = SiteOperator(Basis, QN(-0.5,3), LatticeCommute::Bosonic);
   Um = SiteOperator(Basis, QN(0.5,-3), LatticeCommute::Bosonic);

   P = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   R = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   I = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);
   U = SiteOperator(Basis, QN(0,0), LatticeCommute::Bosonic);

   L3("1", "1")   =  1.0;
   L3("-1", "-1") = -1.0;

   L8("1", "1")   =  1.0 / std::sqrt(3.0);
   L8("-1", "-1") =  1.0 / std::sqrt(3.0);
   L8("0", "0")   = -2.0 / std::sqrt(3.0);

   Sz("1", "1")   =  0.5;
   Sz("-1", "-1") = -0.5;

   Qz("1", "1")   =  1.0;
   Qz("-1", "-1") =  1.0;
   Qz("0", "0")   = -2.0;

   Tp("1", "-1")  =  2.0;
   Tm("-1", "1")  =  2.0;

   Vp("1", "0")   =  2.0;
   Vm("0", "1")   =  2.0;

   Up("-1", "0")  =  2.0;
   Um("0", "-1")  =  2.0;

   I("1", "1")    =  1.0;
   I("-1", "-1")  =  1.0;
   I("0", "0")    =  1.0;

   P = I;
   R = I;

   U("1", "1")    = -1.0;
   U("-1", "-1")  = -1.0;
   U("0", "0")    =  1.0;

   Site["I"] = I;
   Site["P"] = P;
   Site["R"] = R;
   Site["U"] = U;
   Site["L3"] = L3;
   Site["L8"] = L8;
   Site["Sz"] = Sz;
   Site["Qz"] = Qz;
   Site["Tp"] = Tp;
   Site["Tm"] = Tm;
   Site["Vp"] = Vp;
   Site["Vm"] = Vm;
   Site["Up"] = Up;
   Site["Um"] = Um;

   Site.set_operator_descriptions(OpDescriptions);

   return Site;
}
