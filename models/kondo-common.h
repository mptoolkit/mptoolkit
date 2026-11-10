// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/kondo-common.h
//
// Copyright (C) 2026 Ian McCulloch <ian@qusim.net>
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

#ifndef MPTOOLKIT_MODELS_KONDO_COMMON_H
#define MPTOOLKIT_MODELS_KONDO_COMMON_H

#include "lattice/latticesite.h"
#include "quantumnumbers/su2.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

// half_int's stream extraction rounds, so validate the command-line value
// before converting it.  In particular, --Spin 0.7 must not silently mean 1/2.
inline half_int KondoSpin(double Spin)
{
   if (!std::isfinite(Spin) || Spin < 0 ||
       Spin > (std::numeric_limits<int>::max()-2)/2 || std::floor(2*Spin) != 2*Spin)
      throw std::invalid_argument("Spin must be a nonnegative integer or half-integer");
   return half_int(Spin);
}

namespace KondoDetail
{

// Couple the electron first and the local moment second.  This order fixes
// the Clebsch-Gordan phases to agree with the original spin-1/2 Kondo site.
inline LatticeSite MakeSite(LatticeSite const& Electron, half_int Spin,
                           std::string const& SpinSym, std::string const& Description,
                           std::string const& NumberName = "N")
{
   if (Spin < 0)
      throw std::invalid_argument("Spin must be nonnegative");

   SymmetryList Symmetry = Electron.GetSymmetryList();
   QuantumNumbers::QNConstructor<QuantumNumbers::SU2> SQN(SpinSym+":SU(2)");
   SiteBasis MomentBasis(Symmetry);
   MomentBasis.push_back("moment", CoerceSymmetryList(SQN(Spin), Symmetry));
   SiteOperator MomentI(MomentBasis, QuantumNumber(Symmetry), LatticeCommute::Bosonic);
   SiteOperator MomentS(MomentBasis, CoerceSymmetryList(SQN(1), Symmetry), LatticeCommute::Bosonic);
   MomentI(0,0) = 1;
   if (Spin > 0)
      MomentS(0,0) = std::sqrt(Spin.to_double()*(Spin.to_double()+1));

   SiteProductBasis Product(Electron["I"].Basis(), MomentBasis);
   std::vector<int> Order;
   for (unsigned i = 0; i < Product.size(); ++i)
      Order.push_back(i);
   std::sort(Order.begin(), Order.end(), [&](int i, int j) {
      auto a = Product.PBasis().rmap(i).first;
      auto b = Product.PBasis().rmap(j).first;
      return a != b ? a < b : degree(Product.Basis().qn(i)) < degree(Product.Basis().qn(j));
   });

   SiteBasis Basis(Symmetry);
   for (int i : Order)
   {
      std::string Label = Electron["I"].Basis().Label(Product.PBasis().rmap(i).first);
      if (Label == "single")
      {
         bool Lower = degree(Product.Basis().qn(i)) < 2*Spin.to_double()+1;
         Label = Spin == half_int(0.5) ? (Lower ? "singlet" : "triplet")
                                      : (Lower ? "lower" : "upper");
      }
      Basis.push_back(Label, Product.Basis().qn(i));
   }

   auto Lift = [&](SiteOperator const& E, SiteOperator const& F) {
      SiteOperator const Op = tensor_prod(E, F);
      SiteOperator Result(Basis, Op.TransformsAs(), Op.Commute());
      for (unsigned i = 0; i < Order.size(); ++i)
         for (unsigned j = 0; j < Order.size(); ++j)
         {
            std::complex<double> Value = Op(Order[i], Order[j]);
            if (Value != 0.0)
               Result(i,j) = Value;
         }
      return Result;
   };

   LatticeSite Site(Description);
   for (auto Name : {"I", "P", "R", "Hu", "Pdouble", "Pg", "C", "CH", "N_S", "N_H"})
      Site[Name] = Lift(Electron[Name], MomentI);
   Site["N"] = Lift(Electron[NumberName], MomentI);
   Site["Ep"] = Lift(Electron["Qp"], MomentI);
   Site["Em"] = Lift(Electron["Qm"], MomentI);
   // With our spherical-spinor conventions Ep = 2*c_up^dagger*c_down^dagger
   // and Em = -adjoint(Ep).  Keep these legacy operators unchanged, and
   // provide canonically normalized pairs for superconducting Hamiltonians.
   Site["PairH"] = 0.5*Site["Ep"];
   Site["Pair"] = adjoint(Site["PairH"]);
   Site["Sc"] = Lift(Electron[SpinSym], MomentI);
   Site["Sf"] = Lift(Electron["I"], MomentS);
   Site["S"] = Site["Sc"] + Site["Sf"];
   Site["ScSf"] = -std::sqrt(3.0)*prod(Site["Sc"], Site["Sf"], QuantumNumber(Symmetry));
   Site.arg("Spin") = Spin.to_double();

   OperatorDescriptions Descriptions;
   Descriptions.add_operators()
      ("I", "identity")
      ("R", "reflection")
      ("C", "annihilation operator")
      ("CH", "creation operator")
      ("P", "fermion parity")
      ("N", "electron number operator")
      ("N_S", "projector onto single electron occupancy")
      ("N_H", "projector onto empty or double electron occupancy")
      ("S", "total spin vector Sc + Sf")
      ("Sc", "conduction spin vector")
      ("Sf", "local moment spin vector")
      ("Ep", "create double-occupied site (legacy eta-pairing normalization)")
      ("Em", "annihilate double-occupied site (legacy eta-pairing normalization)")
      ("PairH", "onsite singlet pair creation c_up^dagger c_down^dagger")
      ("Pair", "onsite singlet pair annihilation c_down c_up")
      ("Hu", "symmetrized Coulomb operator (n_up - 1/2) * (n_down - 1/2)")
      ("Pdouble", "projector onto double electron occupancy")
      ("ScSf", "Kondo spin interaction Sc . Sf")
      ("Pg", "Gutzwiller projector = 1-Pdouble");
   Site.set_operator_descriptions(Descriptions);
   return Site;
}

} // namespace KondoDetail
#endif
