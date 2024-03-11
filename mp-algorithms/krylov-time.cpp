// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp-algorithms/krylov-time.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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

#include "krylov-time.h"
#include "linearalgebra/eigen.h"
#include "tensor/tensor_eigen.h"
#include "matrixproduct/matrixproduct-sum.h"
#include "common/messagelogger.h"
#include "linearalgebra/matrix_utility.h"

using MessageLogger::Logger;
using MessageLogger::msg_log;

PStream::opstream& operator<<(PStream::opstream& out, KrylovSolver const& s)
{
   return out << s.Psi
              << s.Krylov
              << s.Ident
              << s.Timestep
      ;
}

PStream::ipstream& operator>>(PStream::ipstream& in, KrylovSolver& s)
{
   return in >> s.Psi
             >> s.Krylov
             >> s.Ident
             >> s.Timestep
      ;
}

KrylovSolver::KrylovSolver(MPWavefunction const& Psi_,
                           MPOperator const& H_,
                           std::complex<double> Timestep_,
                           int NumKrylovVectors_)
   : Psi(Psi_), Krylov(H_, std::vector<MPWavefunction>(NumKrylovVectors_+1, Psi_)),
     Ident(Psi_.GetSymmetryList()), Timestep(Timestep_)
{
}

void KrylovSolver::ShiftLeftAndExpand()
{
   Krylov.ShiftLeftAndExpand();
}

void KrylovSolver::ShiftRightAndExpand()
{
   Krylov.ShiftRightAndExpand();
}

void KrylovSolver::ExpandLeft()
{
   Krylov.ExpandLeft();
}

void KrylovSolver::ExpandRight()
{
   Krylov.ExpandRight();
}

void KrylovSolver::UpdateKrylovBasis()
{
   Krylov.ConstructKrylovBasis();
}

void KrylovSolver::TruncateLeft(int MaxStates, double CFactor)
{
   MatrixOperator RhoL = scalar_prod(Krylov.Krylov[0].Center(), herm(Krylov.Krylov[0].Center()));
   MatrixOperator RhoR = scalar_prod(herm(Krylov.Krylov[0].Center()), Krylov.Krylov[0].Center());
   Krylov.TruncateLeft(CFactor, RhoL, RhoR, MaxStates);
}

void KrylovSolver::TruncateRight(int MaxStates, double CFactor)
{
   MatrixOperator RhoL = scalar_prod(Krylov.Krylov[0].Center(), herm(Krylov.Krylov[0].Center()));
   MatrixOperator RhoR = scalar_prod(herm(Krylov.Krylov[0].Center()), Krylov.Krylov[0].Center());
   Krylov.TruncateRight(CFactor, RhoL, RhoR, MaxStates);
}

void KrylovSolver::AdvanceTimestep(int MaxStates)
{
   Matrix<std::complex<double> > M = Krylov.sub_H;
   Vector<std::complex<double> > Eigen = Timestep * DiagonalizeHermitian(M);
   Eigen = transform(Eigen, LinearAlgebra::Exp<std::complex<double> >());
   Matrix<std::complex<double> > X = herm(M) * diagonal_matrix(Eigen) * M;
   Vector<std::complex<double> > v = X(0, LinearAlgebra::all);
   std::vector<MPWavefunction> K = Krylov.Krylov;
   K.pop_back();
   for (std::size_t i = 0; i < K.size(); ++i)
   {
      K[i] *= v[i];
   }

   Psi = fancy_sum(K, MaxStates);
   Krylov.ResetK0(Psi);
}

void KrylovSolver::DebugCheckBasis() const
{
   Krylov.DebugCheckBasis();
}
