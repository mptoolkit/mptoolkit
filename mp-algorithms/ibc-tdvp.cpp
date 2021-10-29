// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/ibc-tdvp.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
// Copyright (C) 2021 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#include "ibc-tdvp.h"
//#include "lanczos-exponential.h"
#include "triangular_mpo_solver.h"
#include "tensor/regularize.h"
#include "tensor/tensor_eigen.h"
#include "linearalgebra/eigen.h"
#include "linearalgebra/exponential.h"

IBC_TDVP::IBC_TDVP(IBCWavefunction const& Psi_, BasicTriangularMPO const& Ham_,
                   std::complex<double> Timestep_, int MaxIter_, double ErrTol_,
                   double GMRESTol_, StatesInfo SInfo_, int Verbose_)
   : GMRESTol(GMRESTol_)
{
   // TODO: Fix member initializer list.
   Hamiltonian = Ham_;
   Timestep = Timestep_;
   MaxIter = MaxIter_;
   ErrTol = ErrTol_;
   SInfo = SInfo_;
   Verbose = Verbose_;

   // Initialize Psi and Ham.
   MatrixOperator Lambda;
   std::tie(Psi, Lambda) = get_left_canonical(Psi_.Window);

   PsiLeft = Psi_.Left;
   PsiRight = Psi_.Right;

   // Set Hamiltonian sizes to match the unit cell/window sizes.
   // TODO: Is there a better place to do this?
   HamiltonianLeft = std::move(Hamiltonian);
   HamiltonianRight = std::move(Hamiltonian);

   if (Hamiltonian.size() < Psi.size())
      Hamiltonian = repeat(Hamiltonian, Psi.size() / Hamiltonian.size());

   if (HamiltonianLeft.size() < PsiLeft.size())
      HamiltonianLeft = repeat(HamiltonianLeft, PsiLeft.size() / HamiltonianLeft.size());

   if (HamiltonianRight.size() < PsiRight.size())
      HamiltonianRight = repeat(HamiltonianRight, PsiRight.size() / HamiltonianRight.size());

   if (Verbose > 0)
      std::cout << "Constructing Hamiltonian block operators..." << std::endl;

   StateComponent BlockHamL = Initial_E(HamiltonianLeft, PsiLeft.Basis2());
   std::complex<double> LeftEnergy = SolveSimpleMPO_Left(BlockHamL, PsiLeft, HamiltonianLeft,
                                                         GMRESTol, Verbose-1);
   if (Verbose > 0)
      std::cout << "Starting energy (left eigenvalue) = " << LeftEnergy << std::endl;

   StateComponent BlockHamR = Initial_F(HamiltonianRight, PsiRight.Basis1());
   std::complex<double> RightEnergy = SolveSimpleMPO_Right(BlockHamR, PsiRight, HamiltonianRight,
                                                           GMRESTol, Verbose-1);
   if (Verbose > 0)
      std::cout << "Starting energy (right eigenvalue) = " << RightEnergy << std::endl;

   Ham.push_left(BlockHamL);
   Ham.push_right(BlockHamR);

   C = Psi.begin();
   H = Hamiltonian.begin();
   while (C != Psi.end())
   {
      if (Verbose > 1)
         std::cout << "Site " << (Ham.size_left()) << std::endl;
      Ham.push_left(contract_from_left(*H, herm(*C), Ham.left(), *C));
      MaxStates = std::max(MaxStates, (*C).Basis2().total_dimension());
      ++H, ++C;
   }

   Psi.set_back(prod(Psi.get_back(), Lambda));

   // Initialize to the right-most site.
   Ham.pop_left();
   C = Psi.end();
   --C;
   --H;
   Site = Psi.size() - 1;

   LeftStop = 0;
   RightStop = Psi.size() - 1;
}

IBCWavefunction
IBC_TDVP::Wavefunction() const
{
   MatrixOperator I = MatrixOperator::make_identity(Psi.Basis2());
   WavefunctionSectionLeft PsiWindow = WavefunctionSectionLeft::ConstructFromLeftOrthogonal(std::move(Psi), I, Verbose);

   return IBCWavefunction(PsiLeft, PsiWindow, PsiRight);
}
