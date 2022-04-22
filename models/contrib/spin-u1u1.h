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

// SU(3) spin chain in the fundamental representation,
// broken down to U(1)xU(1)
//
// This uses the basis from https://arxiv.org/abs/0801.3565
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
// We take the U(1) generators to be L3/2 and L8 * sqrt(3)/2.  We denote the basis by the eigenvalues of L3,
// being (1,-1,0), which have quantum numbers (1/2,1/2), (-1/2,1/2) and (0,-1) respectively.
//
// A conventient representation is
// Tp = L1 + i*L2
// Tm = L1 - i*L2
// Vp = L4 + i*L5
// Vp = L4 - i*L5
// Up = L6 + i*L7
// Um = L6 - i*L7
//
// Tp,Tm are raising/lowering operators in the first two states.
// Up,Um are raising/lowering operators in the second two states.
// Vp,Vm are raising/lowering operators between the first and last states.
//
// Relation so S=1 matrices
// ========================
//
// We construct the S=1 Pauli matrices, using the same basis of (1,-1,0).
//
// Sz = [  1  0  0 ]
//      [  0 -1  0 ]
//      [  0  0  0 ]
//
// Sx = (1/sqrt(2)) [  0  0  1 ]
//                  [  0  0  1 ]
//                  [  1  1  0 ]
//
// Sy = (1/sqrt(2)) [  0  0 -i ]
//                  [  0  0  i ]
//                  [  i -i  0 ]
//
// Sp = sqrt(2) [  0  0  1 ]
//              [  0  0  0 ]
//              [  0  1  0 ]
//
// Sm = sqrt(2) [  0  0  0 ]
//              [  0  0  1 ]
//              [  1  0  0 ]
//
// Sz = L3
// Sx = (1/sqrt(2)) (L4 + L6)
// Sy = (1/sqrt(2)) (L5 - L7)
//
// Hamiltonians
// ============
//
// The SU(3) spin chain has the Hamiltonian
//
// H_{SU(3)} = (J/4) sum_j sum_{a=1}^8 La(j) La(j+1)
//
// This is closely related to the S=1 bilinear biquadratic model
// H_{bq} = J sum_j S_j \cdot S_{j+1} + (S_j \cdot S_{j+1})^2
//
// via S_j \cdot S_{j+1} + (S_j \cdot S_{j+1})^2 = 1 + 1/3 + (1/2) sum_{a=1}^8 La(j) La(j+1)
//
// Hence the translation is
// H_{bq} = (4/3) N + 2 H_{SU(3)}
//
// The energy per site of the bilinear biquadriatic model (with J=1) is known to be
// E_BQ = -ln(3) - pi/(3 sqrt(3)) + 2 = 0.296787923253818
//
// so the energy per site of the SU(3) spin chain (with J=1) is
// E_SU3 = 1/3 - ln(3)/2 - pi/(6 sqrt(3)) = -0.518272705039758
//
// To write the Hamiltonian in terms of the generators, we can use
//
// H_{bq} = (4/3) + (1/2) [ (1/2) (Tp Tm + Vp Vm + Up Um) + L3 L3 + L8 L8 ]
//
// and H_{SU(3)} = (1/4) [ (1/2) (Tp Tm + Vp Vm + Up Um) + L3 L3 + L8 L8 ]
//

inline
LatticeSite SpinU1U1()
{
   SymmetryList Symmetry("Sz:U(1),Qz:U(1)");
   QuantumNumbers::QNConstructor<QuantumNumbers::U1, QuantumNumbers::U1>
      QN(Symmetry);
   SiteBasis Basis(Symmetry);
   SiteOperator L3, L8;
   SiteOperator Tp, Tm, Up, Um, Vp, Vm, Sz, Qz;

   SiteOperator R, P, I, U;
   LatticeSite Site;

   // Basis is labelled by Sz = L3/2 and Qz = L8 * sqrt(3)/2
   // we label the local basis as 1,-1,0, the eigenvalue of L3
   Basis.push_back("1", QN(0.5, 0.5));
   Basis.push_back("-1", QN(-0.5, 0.5));
   Basis.push_back("0", QN(0, -1));

   OperatorDescriptions OpDescriptions;
   OpDescriptions.add_operators()
      ("I"   , "identity")
      ("R"   , "reflection")
      ("P"   , "fermion parity")
      ("U"   , "unitary rotation in lambda_8 direction")
      ("L3"  , "U(1) generator lambda_3")
      ("L8"  , "U(1) generator lambda_8")
      ("Sz"  , "U(1) generator lambda_3 scaled to have integer eigenvalues -1/2,1/2")
      ("Qz"  , "U(1) generator scaled to have integer eigenvalues 1/2,1/2,-1")
      ("Tp"  , "T raising operator")
      ("Tm"  , "T lowering operator")
      ("Vp"  , "V raising operator")
      ("Vm"  , "V lowering operator")
      ("Up"  , "U raising operator")
      ("Um"  , "U lowering operator")
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
   Vp = SiteOperator(Basis, QN(0.5,1.5), LatticeCommute::Bosonic);
   Vm = SiteOperator(Basis, QN(-0.5,-1.5), LatticeCommute::Bosonic);

   // U+ = L6 + i*L7
   // U- = L6 - i*L7
   Up = SiteOperator(Basis, QN(-0.5,1.5), LatticeCommute::Bosonic);
   Um = SiteOperator(Basis, QN(0.5,-1.5), LatticeCommute::Bosonic);

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

   Qz("1", "1")   =  0.5;
   Qz("-1", "-1") =  0.5;
   Qz("0", "0")   = -1.0;

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
