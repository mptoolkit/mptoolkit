// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// tests/helpers/test-kondo.cpp
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
#include "models/kondo-u1su2.h"
#include "models/kondo-su2.h"
#include "linearalgebra/eigen.h"
#include <iostream>

using Dense = LinearAlgebra::Matrix<std::complex<double>>;

void Require(bool Pass, std::string const& What)
{
   if (!Pass)
      throw std::runtime_error(What);
}

void Close(Dense const& A, Dense const& B, std::string const& What)
{
   Require(size1(A) == size1(B) && size2(A) == size2(B), What+": dimensions");
   Require(norm_frob(A-B) < 2E-11, What+": matrix mismatch");
}

Dense Reduced(SiteOperator const& Op)
{
   Dense M(Op.Basis().size(), Op.Basis().size(), 0.0);
   for (unsigned i = 0; i < size1(M); ++i)
      for (unsigned j = 0; j < size2(M); ++j)
         M(i,j) = Op(i,j);
   return M;
}

struct State
{
   int Multiplet;
   QuantumNumber Q;
   Projection M;
};

std::vector<State> States(SiteBasis const& Basis)
{
   std::vector<State> Result;
   for (unsigned i = 0; i < Basis.size(); ++i)
   {
      std::vector<Projection> M;
      enumerate_projections(Basis.qn(i), std::back_inserter(M));
      for (auto const& m : M)
         Result.push_back({int(i), Basis.qn(i), m});
   }
   return Result;
}

// Wigner-Eckart expansion, independent of the tensor-product/9j machinery
// used to build the model.  This helper is for the SU(2)-only site.
Dense Expand(SiteOperator const& Op, double Component = 0)
{
   auto StatesList = States(Op.Basis());
   Projection K(Op.GetSymmetryList(), ConvertToString(Component));
   Dense M(StatesList.size(), StatesList.size(), 0.0);
   for (unsigned i = 0; i < StatesList.size(); ++i)
      for (unsigned j = 0; j < StatesList.size(); ++j)
      {
         auto const& a = StatesList[i];
         auto const& b = StatesList[j];
         if (is_transform_target(b.Q, Op.TransformsAs(), a.Q))
            M(i,j) = Op(a.Multiplet,b.Multiplet)*
               clebsch_gordan(b.Q, Op.TransformsAs(), a.Q, b.M, K, a.M);
      }
   return M;
}

void CheckSpin(double s)
{
   auto Site = KondoSU2(half_int(s));
   auto ChargeSite = KondoU1SU2(half_int(s));
   auto Custom = KondoU1SU2(half_int(s), "charge", "spin");
   auto CustomSU2 = KondoSU2(half_int(s), "spin");
   auto const& Basis = Site["I"].Basis();
   auto Coupled = States(Basis);
   int d = int(2*s+1), n = 4*d;
   Require(Coupled.size() == unsigned(n), "physical site dimension");
   Require(Basis.size() == (s == 0 ? 3 : 4), "multiplet count");

   for (auto It = Site.begin_operator(); It != Site.end_operator(); ++It)
   {
      auto const& Name = It->first;
      Close(Reduced(It->second), Reduced(ChargeSite[Name]), Name+": charge symmetry");
      Close(Reduced(It->second), Reduced(Custom[Name]), Name+": custom symmetry names");
      Close(Reduced(It->second), Reduced(CustomSU2[Name]), Name+": custom SU2 name");
   }
   Require(Site["C"].Commute() == LatticeCommute::Fermionic, "C statistics");
   Require(Site["CH"].Commute() == LatticeCommute::Fermionic, "CH statistics");
   Require(Site["Pair"].Commute() == LatticeCommute::Bosonic, "pair statistics");

   // Explicit uncoupled basis: |empty,up,down,double> x |m=-S,...,S>.
   Dense U(n,n,0.0), I(n,n,0.0), Cu(n,n,0.0), Cd(n,n,0.0);
   Dense Fz(n,n,0.0), Fp(n,n,0.0);
   SymmetryList SL("S:SU(2)");
   QuantumNumbers::QNConstructor<QuantumNumbers::SU2> QN(SL);
   for (int f = 0; f < d; ++f)
   {
      Cu(f,d+f) = 1;
      Cu(2*d+f,3*d+f) = 1;
      Cd(f,2*d+f) = 1;
      Cd(d+f,3*d+f) = -1;
      double m = f-s;
      for (int e = 0; e < 4; ++e)
      {
         I(e*d+f,e*d+f) = 1;
         Fz(e*d+f,e*d+f) = m;
         if (f+1 < d)
            Fp(e*d+f+1,e*d+f) = std::sqrt((s-m)*(s+m+1));
         for (int k = 0; k < n; ++k)
         {
            std::string Label = Basis.Label(Coupled[k].Multiplet);
            if ((e == 0 && Label == "empty") || (e == 3 && Label == "double"))
               U(e*d+f,k) = Coupled[k].M == Projection(SL,ConvertToString(m)) ? 1.0 : 0.0;
            if ((e == 1 || e == 2) && Label != "empty" && Label != "double")
               U(e*d+f,k) = clebsch_gordan(QN(0.5), QN(s), Coupled[k].Q,
                  Projection(SL, e == 1 ? "0.5" : "-0.5"),
                  Projection(SL,ConvertToString(m)), Coupled[k].M);
         }
      }
   }
   Close(Dense(herm(U)*U), I, "coupled basis unitarity");
   Dense Nu = herm(Cu)*Cu, Nd = herm(Cd)*Cd;
   Dense N = Nu+Nd, D = Nu*Nd, P = I-2.0*N+4.0*D;
   Dense Cz = 0.5*(Nu-Nd), Cp = herm(Cu)*Cd;
   Dense K = Cz*Fz+0.5*(Cp*herm(Fp)+herm(Cp)*Fp);
   Dense PairH = herm(Cu)*herm(Cd), Pair = Cd*Cu;
   auto Check = [&](std::string const& Name, Dense const& Expected, double Component = 0) {
      Close(Expand(Site[Name],Component), Dense(herm(U)*Expected*U), Name);
   };
   Check("I", I);
   Check("C", Cu, -0.5);
   Check("C", Dense(-Cd), 0.5);
   Check("CH", Dense(herm(Cu)), 0.5);
   Check("CH", Dense(herm(Cd)), -0.5);
   Check("N", N);
   Check("P", P);
   Check("Pdouble", D);
   Check("R", Dense(I-2.0*D));
   Check("Hu", Dense(D-0.5*N+0.25*I));
   Check("Pg", Dense(I-D));
   Check("N_S", Dense(N-2.0*D));
   Check("N_H", Dense(I-N+2.0*D));
   Check("Sc", Cz);
   Check("Sc", Dense((-1/std::sqrt(2.0))*Cp), 1);
   Check("Sf", Fz);
   Check("Sf", Dense((-1/std::sqrt(2.0))*Fp), 1);
   Check("S", Dense(Cz+Fz));
   Check("S", Dense((-1/std::sqrt(2.0))*(Cp+Fp)), 1);
   Check("ScSf", K);
   Check("Ep", Dense(2.0*PairH));
   Check("Em", Dense(-2.0*Pair));
   Check("PairH", PairH);
   Check("Pair", Pair);

   // The central SU(2) element fixes electron parity, including the offset
   // from half-integer local moments.  Empty/double are NOT spinless for S>0.
   for (unsigned i = 0; i < Basis.size(); ++i)
   {
      int TwiceJ = degree(Basis.qn(i))-1;
      double Parity = (TwiceJ-int(2*s))%2 == 0 ? 1 : -1;
      Require(std::abs(std::complex<double>(Site["P"](i,i))-Parity) < 1E-12, "SU2 parity offset");
   }

   // Analytic spectrum with both real and imaginary pairing and nonzero mu.
   double mu = 0.37, JK = 0.83;
   std::complex<double> Delta(0.41,-0.29);
   Dense H = JK*Expand(Site["ScSf"])-mu*Expand(Site["N"])
      +Delta*Expand(Site["PairH"])+std::conj(Delta)*Expand(Site["Pair"]);
   Close(H, Dense(herm(H)), "paired Hamiltonian Hermiticity");
   auto Eigenvalues = EigenvaluesHermitian(H);
   std::vector<double> Expected;
   for (int i = 0; i < d; ++i)
   {
      Expected.push_back(-mu-std::sqrt(mu*mu+std::norm(Delta)));
      Expected.push_back(-mu+std::sqrt(mu*mu+std::norm(Delta)));
   }
   for (int i = 0; i < int(2*s); ++i)
      Expected.push_back(-mu-JK*(s+1)/2);
   for (int i = 0; i < int(2*s+2); ++i)
      Expected.push_back(-mu+JK*s/2);
   std::sort(Expected.begin(),Expected.end());
   for (int i = 0; i < n; ++i)
      Require(std::abs(Eigenvalues[i]-Expected[i]) < 2E-11, "paired Hamiltonian spectrum");
}

void CheckLegacy()
{
   auto Site = KondoU1SU2();
   std::vector<std::string> Labels{"empty","double","singlet","triplet"};
   for (int i = 0; i < 4; ++i)
      Require(Site["I"].Basis().Label(i) == Labels[i], "legacy basis labels/order");
   Dense C(4,4,0.0), Sc(4,4,0.0), Sf(4,4,0.0), S(4,4,0.0);
   C(0,2) = std::sqrt(0.5); C(0,3) = std::sqrt(1.5); C(2,1) = -1; C(3,1) = 1;
   Sc(2,3) = -std::sqrt(0.75); Sc(3,2) = 0.5; Sc(3,3) = std::sqrt(0.5);
   Sf(0,0) = Sf(1,1) = std::sqrt(0.75);
   Sf(2,3) = std::sqrt(0.75); Sf(3,2) = -0.5; Sf(3,3) = std::sqrt(0.5);
   S(0,0) = S(1,1) = std::sqrt(0.75); S(3,3) = std::sqrt(2.0);
   Close(Reduced(Site["C"]), C, "legacy C");
   Close(Reduced(Site["Sc"]), Sc, "legacy Sc");
   Close(Reduced(Site["Sf"]), Sf, "legacy Sf");
   Close(Reduced(Site["S"]), S, "legacy total S");
   auto Custom = KondoU1SU2("charge", "spin");
   Close(Reduced(Custom["C"]), C, "legacy constructor overload");
}

int main()
{
   try
   {
      CheckLegacy();
      for (double s : {0.0, 0.5, 1.0, 1.5, 2.0, 3.5})
         CheckSpin(s);
      for (double s : {-0.5, 0.7, std::numeric_limits<double>::infinity(),
                       std::numeric_limits<double>::quiet_NaN(), 1E30})
      {
         bool Rejected = false;
         try { KondoSpin(s); }
         catch (std::invalid_argument const&) { Rejected = true; }
         Require(Rejected, "invalid spin accepted");
      }
      std::cout << "Kondo site matrices, legacy conventions, parity and paired spectra passed\n";
   }
   catch (std::exception const& e)
   {
      std::cerr << e.what() << '\n';
      return 1;
   }
}
