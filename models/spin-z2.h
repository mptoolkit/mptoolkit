// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/spin-z2.h
//
// Copyright (C) 2004-2018 Ian McCulloch <ian@qusim.net>
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
#include "quantumnumbers/z2.h"




std::string StateName(half_int s, bool Symmetric)
{
   std::string Result = boost::lexical_cast<std::string>(s);
   Result += Symmetric ? 's' : 'a';
   return Result;
}

inline
LatticeSite SpinZ2(half_int Spin)
{
   SymmetryList Symmetry("Z:Z2");
   QuantumNumbers::QNConstructor<QuantumNumbers::Z2> QN(Symmetry);
   SiteBasis Basis(Symmetry);
   SiteOperator Sz, Sx, Sy, mSz, R, P, I, Z;
   LatticeSite Site("Spin "+to_string_fraction(Spin) + " " + Symmetry.FullName());

   PRECONDITION(Spin >= 0)("Spin must be >= 0")(Spin);

   half_int Base = 0.5; // lowest possible quantum number
   if (Spin.is_integral())
   {
      Base = 1;
      Basis.push_back(StateName(0,true), QN(1));
   }
   for (half_int s = Base; s <= Spin; ++s)
   {
      Basis.push_back(StateName(s, true), QN(1));
      Basis.push_back(StateName(s, false), QN(-1));
   }

   Sx = SiteOperator(Basis, QN(1), LatticeCommute::Bosonic);
   Sy = SiteOperator(Basis, QN(-1), LatticeCommute::Bosonic);
   Sz = SiteOperator(Basis, QN(-1), LatticeCommute::Bosonic);
   P = SiteOperator(Basis, QN(1), LatticeCommute::Bosonic);
   R = SiteOperator(Basis, QN(1), LatticeCommute::Bosonic);
   I = SiteOperator::Identity(Basis);

   P = I; // bosonic, so the parity operator is identity
   R = I; // spatial reflection


   // Set the matrix elements for Sz
   // We only set half the matrix elements of Sx,Sy, and at the end
   // make them hermitian
   // Since Base > 0 we do the s=0 components later (these are special anyway)
   for (half_int s = Base; s <= Spin; ++s)
   {
      Sz(StateName(s, true), StateName(s, false)) = s.to_double();
      Sz(StateName(s, false), StateName(s, true)) = s.to_double();

      if (s < Spin)
      {
         Sx(StateName(s, false), StateName(s+1, false))
            = 0.5 * sqrt((Spin-s)*(Spin+s+1));
         Sx(StateName(s, true), StateName(s+1, true))
            = 0.5 * sqrt((Spin-s)*(Spin+s+1));

         Sy(StateName(s, false), StateName(s+1, true))
            = std::complex<double>(0.0, -0.5 * sqrt((Spin-s)*(Spin+s+1)));
         Sy(StateName(s, true), StateName(s+1, false))
            = std::complex<double>(0.0, -0.5 * sqrt((Spin-s)*(Spin+s+1)));
      }
   }

   // do the s == 0 part
   if (Spin.is_integral())
   {
      Sz("0s", "0s") = 1;
      if (Spin > 0)
      {
         Sx("0s", "1s")
            = sqrt(Spin*(Spin+1)/2);

         Sy("0s", "1a")
            = std::complex<double>(0.0, -sqrt(Spin*(Spin+1)/2));
      }
   }

   Sx = Sx + adjoint(Sx);
   Sy = Sy + adjoint(Sy);

   if (!Spin.is_integral())
   {
      // Some of the matrix elements of Sx,Sy for s=0.5 are diagonal in m
      Sx(StateName(0.5, true), StateName(0.5, true))
         = 0.5 * (Spin+0.5);
      Sx(StateName(0.5, false), StateName(0.5, false))
         = -0.5 * (Spin+0.5);

      Sy(StateName(0.5, true), StateName(0.5, false))
         = std::complex<double>(0.0, -0.5 * (Spin+0.5));
      Sy(StateName(0.5, false), StateName(0.5, true))
         = std::complex<double>(0.0, 0.5 * (Spin+0.5));
   }

#if 0
   Sx("s", "s") =  0.5;
   Sx("a", "a") = -0.5;

   Sz("s", "a") =  0.5;
   Sz("a", "s") =  0.5;

   Sy("s", "a") = std::complex<double>(0.0, -0.5);
   Sy("a", "s") = std::complex<double>(0.0,  0.5);
#endif

   Site["I"] = I;
   Site["P"] = P;
   Site["R"] = R;
   Site["Sx"] = Sx;
   Site["Sy"] = Sy;
   Site["Sz"] = Sz;
   Site["X"] = 2*Sx;
   Site["Y"] = 2*Sy;
   Site["Z"] = 2*Sz;
   Site["Sz2"] = Sz * Sz;
   return Site;
}
