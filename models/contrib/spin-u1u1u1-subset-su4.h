// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/contrib/spin-u1u1u1-subset-su4.h
//
// Copyright (C) 2004-2022 Ian McCulloch <ian@qusim.net>
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
// In terms of the basis states |0> |1> |2> |3>, these are:
//
// Tp|0> = 2 |1>
// Vp|0> = 2 |2>
// Up|1> = 2 |2>
// Wp|0> = 2 |3>
// Xp|1> = 2 |3>
// Yp|2> = 2 |3>
//
// We can assemble the raising and lowering operators into SU(2) generators by adding:
// Tz = L3
// Vz = ( sqrt(3)*L8 + L3 ) / 2
// Uz = ( sqrt(3)*L8 - L3 ) / 2
// Wz = ( sqrt(6)*L15 + (sqrt(3)/2)*L8 + (3/2)*L3 ) / 3
// Xz = ( sqrt(6)*L15 + (sqrt(3)/2)*L8 - (3/2)*L3 ) / 3
// Yz = ( sqrt(6)*L15 - sqrt(3)*L8 ) / 3
//
// These generators differ from the usual spin-1/2 generators as a factor 2 larger.
//
// These split into SU(2)xSU(2) in 3 different ways, with
// [T,Y] = 0
// [V,X] = 0
// [U,W] = 0
//
// The SU(4) Hamiltonian is
// H_{SU(4)} = (1/4) [ (1/2) (Tp Tm + Vp Vm + Up Um + Wp Wm + Xp Xm + Yp Ym) + L3 L3 + L8 L8 + L15 L15 ]
//

inline
LatticeSite SpinU1U1U1()
{
   SymmetryList Symmetry("L3:U(1),L8:U(1),L15:U(1)");
   QuantumNumbers::QNConstructor<QuantumNumbers::U1, QuantumNumbers::U1, QuantumNumbers::U1>
      QN(Symmetry);
   SiteBasis Basis(Symmetry);
   SiteOperator L3, L8, L15;
   SiteOperator M3, M8, M15;
   SiteOperator Tp, Tm, Up, Um, Vp, Vm, Wp, Wm, Xp, Xm, Yp, Ym;
   SiteOperator Tz, Vz, Uz, Wz, Xz, Yz;

   SiteOperator R, P, I;
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
      ("L3"  , "U(1) generator lambda_3")
      ("L8"  , "U(1) generator lambda_8")
      ("L15" , "U(1) generator lambda_15")
      ("M3"  , "U(1) generator lambda_3 scaled to have integer eigenvalues 1,-1,0,0")
      ("M8"  , "U(1) generator lambda_8 scaled to have integer eigenvalues 1,1,-2,0")
      ("M15" , "U(1) generator lambda_15 scaled to have integer eigenvalues 1,1,1,-3")
      ("Tp"  , "T raising operator")
      ("Tm"  , "T lowering operator")
      ("Tz"  , "T z-spin operator")
      ("Vp"  , "V raising operator")
      ("Vm"  , "V lowering operator")
      ("Vz"  , "T z-spin operator")
      ("Up"  , "U raising operator")
      ("Um"  , "U lowering operator")
      ("Uz"  , "T z-spin operator")
      ("Wp"  , "W raising operator")
      ("Wm"  , "W lowering operator")
      ("Wz"  , "T z-spin operator")
      ("Xp"  , "X raising operator")
      ("Xm"  , "X lowering operator")
      ("Xz"  , "T z-spin operator")
      ("Yp"  , "Y raising operator")
      ("Ym"  , "Y lowering operator")
      ("Yz"  , "T z-spin operator")
      ;

   I = SiteOperator(Basis, QN(0,0,0), LatticeCommute::Bosonic);
   I("0", "0")  =  1.0;
   I("1", "1")  =  1.0;
   I("2", "2")  =  1.0;
   I("3", "3")  =  1.0;
   P = I;
   R = I;

   M3  = SiteOperator(Basis, QN(0,0,0), LatticeCommute::Bosonic);
   M3("0", "0") =  1.0;
   M3("1", "1") = -1.0;
   L3 = M3;

   M8  = SiteOperator(Basis, QN(0,0,0), LatticeCommute::Bosonic);
   M8("0", "0") =  1.0;
   M8("1", "1") =  1.0;
   M8("2", "2") = -2.0;
   L8 = M8 * (1.0 / std::sqrt(3.0));

   M15 = SiteOperator(Basis, QN(0,0,0), LatticeCommute::Bosonic);
   M15("0", "0") =  1.0;
   M15("1", "1") =  1.0;
   M15("2", "2") =  1.0 ;
   M15("3", "3") = -3.0;
   L15 = M15 * (1.0 / std::sqrt(6.0));

   // T+ = L1 + i*L2
   // T- = L1 - i*L2
   // |1><0|
   Tp = SiteOperator(Basis, QN(2,0,0), LatticeCommute::Bosonic);
   Tm = SiteOperator(Basis, QN(-2,0,0), LatticeCommute::Bosonic);
   Tp("0", "1")  =  2.0;
   Tm("1", "0")  =  2.0;

   // V+ = L4 + i*L5
   // V- = L4 - i*L5
   // |2><0|
   Vp = SiteOperator(Basis, QN(1,3,0), LatticeCommute::Bosonic);
   Vm = SiteOperator(Basis, QN(-1,-3,0), LatticeCommute::Bosonic);
   Vp("0", "2")   =  2.0;
   Vm("2", "0")   =  2.0;

   // U+ = L6 + i*L7
   // U- = L6 - i*L7
   // |2><1|
   Up = SiteOperator(Basis, QN(-1,3,0), LatticeCommute::Bosonic);
   Um = SiteOperator(Basis, QN(1,-3,0), LatticeCommute::Bosonic);
   Up("1", "2")  =  2.0;
   Um("2", "1")  =  2.0;

   // W+ = L9 + i*L10
   // W- = L10 - i*L10
   // |3><0|
   Wp = SiteOperator(Basis, QN(1,1,4), LatticeCommute::Bosonic);
   Wm = SiteOperator(Basis, QN(-1,-1,-4), LatticeCommute::Bosonic);
   Wp("0", "3")  =  2.0;
   Wm("3", "0")  =  2.0;

   // X+ = L11 + i*L12
   // X- = L11 - i*L12
   // |3><1|
   Xp = SiteOperator(Basis, QN(-1,1,4), LatticeCommute::Bosonic);
   Xm = SiteOperator(Basis, QN(1,-1,-4), LatticeCommute::Bosonic);
   Xp("1", "3")  =  2.0;
   Xm("3", "1")  =  2.0;

   // Y+ = L13 + i*L14
   // Y- = L13 - i*L14
   // |3><2|
   Yp = SiteOperator(Basis, QN(0,-2,4), LatticeCommute::Bosonic);
   Ym = SiteOperator(Basis, QN(0,2,-4), LatticeCommute::Bosonic);
   Yp("2", "3")  =  2.0;
   Ym("3", "2")  =  2.0;

   Tz = L3;
   Vz = 0.5 * (M8 + L3);
   Uz = 0.5 * (M8 - L3);
   Wz = (1.0/3.0) * (M15 + 0.5*M8 + (3/2)*L3);
   Xz = (1.0/3.0) * (M15 + 0.5*M8 - (3/2)*L3);
   Yz = (1.0/3.0) * (M15 - M8);

   Site["I"] = I;
   Site["P"] = P;
   Site["R"] = R;
   Site["L3"] = L3;
   Site["L8"] = L8;
   Site["L15"] = L15;
   Site["M3"] = M3;
   Site["M8"] = M8;
   Site["M15"] = M15;
   Site["Tp"] = Tp;
   Site["Tm"] = Tm;
   Site["Tz"] = Tz;
   Site["Vp"] = Vp;
   Site["Vm"] = Vm;
   Site["Vz"] = Vz;
   Site["Up"] = Up;
   Site["Um"] = Um;
   Site["Uz"] = Uz;
   Site["Wp"] = Wp;
   Site["Wm"] = Wm;
   Site["Wz"] = Wz;
   Site["Xp"] = Xp;
   Site["Xm"] = Xm;
   Site["Xz"] = Xz;
   Site["Yp"] = Yp;
   Site["Ym"] = Ym;
   Site["Yz"] = Yz;

   Site.set_operator_descriptions(OpDescriptions);

   return Site;
}
