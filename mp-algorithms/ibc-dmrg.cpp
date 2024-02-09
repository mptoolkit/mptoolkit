// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/ibc-dmrg.cpp
//
// Copyright (C) 2024 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#include "ibc-dmrg.h"
#include "mp-algorithms/lanczos.h"

struct HEff1
{
   HEff1(StateComponent const& E_, OperatorComponent const& H_, StateComponent const& F_)
      : E(E_), H(H_), F(F_)
   {
   }

   StateComponent operator()(StateComponent const& x) const
   {
      StateComponent R = operator_prod_inner(H, E, x, herm(F));
      return R;
   }

   StateComponent const& E;
   OperatorComponent const& H;
   StateComponent const& F;
};

struct HEff2
{
   HEff2(StateComponent const& E_, StateComponent const& F_)
      : E(E_), F(F_)
   {
   }

   MatrixOperator operator()(MatrixOperator const& x) const
   {
      MatrixOperator R = operator_prod(E, x, herm(F));
      return R;
   }

   StateComponent const& E;
   StateComponent const& F;
};


IBC_DMRG::IBC_DMRG(IBCWavefunction const& Psi_, WindowHamiltonian const& Ham_, IBC_DMRGSettings const& Settings_)
   : IBC_TDVP(Psi_, Ham_, Settings_)
{
}

void
IBC_DMRG::EvolveCurrentSite(std::complex<double> Tau)
{
   // Evolve current site.
   int Iter = MaxIter;
   double Err = ErrTol;

   double Energy = Lanczos(*C, HEff1(HamL.back(), *H, HamR.front()), Iter, Err);

   if (Verbose > 1)
   {
      std::cout << "Timestep=" << TStep
                << " Site=" << Site
                << " Energy=" << Energy
                << " C Iter=" << Iter
                << " Err=" << Err
                << std::endl;
   }
}

void
IBC_DMRG::IterateLeft(std::complex<double> Tau)
{
   // Perform SVD to right-orthogonalize current site.
   MatrixOperator U;
   RealDiagonalOperator D;

   std::tie(U, D) = OrthogonalizeBasis1(*C);

   // Update the effective Hamiltonian.
   HamR.push_front(contract_from_right(herm(*H), *C, HamR.front(), herm(*C)));

   // Move to the next site.
   --Site;
   --H;
   --C;

   *C = prod(*C, U*D);

   HamL.pop_back();
}

void
IBC_DMRG::IterateRight(std::complex<double> Tau)
{
   // Perform SVD to left-orthogonalize current site.
   MatrixOperator Vh;
   RealDiagonalOperator D;

   std::tie(D, Vh) = OrthogonalizeBasis2(*C);

   // Update the effective Hamiltonian.
   HamL.push_back(contract_from_left(*H, herm(*C), HamL.back(), *C));

   // Move to the next site.
   ++Site;
   ++H;
   ++C;

   *C = prod(D*Vh, *C);

   HamR.pop_front();
}
