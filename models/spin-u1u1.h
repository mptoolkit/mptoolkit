// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/spin-u1u1.h
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
// SU(3) spin chain in the fundamental representation,
// broken down to U(1)xU(1)

#include "lattice/latticesite.h"
#include "quantumnumbers/u1.h"




inline
LatticeSite CreateU1U1SpinSite()
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
   return Site;
}
