// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/split-krylov.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "split-krylov.h"
#include "matrixproduct/wavefunc-utils.h"

PStream::opstream& operator<<(PStream::opstream& out, SplitKrylov const& s)
{
   return out << s.Krylov
              << s.H
              << s.kn1_H_kn
              << s.ki_kj
              << s.sub_H
              << s.Ident
      ;
}

PStream::ipstream& operator>>(PStream::ipstream& in, SplitKrylov& s)
{
   return in >> s.Krylov
             >> s.H
             >> s.kn1_H_kn
             >> s.ki_kj
             >> s.sub_H
             >> s.Ident
      ;
}

SplitKrylov::SplitKrylov(MPOperator const& H_, std::vector<MPWavefunction> const& Krylov_)
   : Krylov(Krylov_), H(H_), kn1_H_kn(Krylov_.size()-1), ki_kj(Krylov_.size(), Krylov_.size()),
     Ident(H_.GetSymmetryList())
{
   this->ResetK0(Krylov[0]);
}

void SplitKrylov::ResetK0(MPWavefunction const& Psi)
{
   // Set K0 and make sure the center matrix is in the same location as the old K0
   int kl = Krylov[0].LeftSize();
   Krylov[0] = Psi;
   Krylov[0].normalize();
   while (Krylov[0].LeftSize() > kl)
      Krylov[0].RotateLeft();
   while (Krylov[0].LeftSize() < kl)
      Krylov[0].RotateRight();

   // Initialize the operator matrix elements for < k(i) | k(j) > and < k(i+1) | H | k(i) >
   for (std::size_t i = 0; i < kn1_H_kn.size(); ++i)
   {
      InitializeSuperblockStack(kn1_H_kn[i], Krylov[i+1], H, Krylov[i]);

      for (std::size_t j = i+1; j < Krylov.size(); ++j)
      {
         InitializeTransformStack(ki_kj(i,j), Krylov[i], Krylov[j]);
      }
   }
}

void SplitKrylov::ShiftRightAndExpand()
{
   // Wavefunctions
   H.RotateRight();
   for (std::size_t i = 0; i < Krylov.size(); ++i)
   {
      Krylov[i].RotateRightExpand();
   }
   // matrix elements of < Krylov[i+1] | H | Krylov[i] >
   for (std::size_t i = 0; i < kn1_H_kn.size(); ++i)
   {
      SuperblockStackRotateRight(kn1_H_kn[i], Krylov[i+1], H, Krylov[i]);
   }
   // matrix elements of <Krylov[i] | Krylov[j] >
   for (std::size_t i = 0; i < Krylov.size()-1; ++i)
   {
      for (std::size_t j = i+1; j < Krylov.size(); ++j)
      {
         TransformStackRotateRight(ki_kj(i,j), Krylov[i], Krylov[j]);
      }
   }
   this->DebugCheckBasis();
}

void SplitKrylov::ShiftLeftAndExpand()
{
   // Ham operator
   H.RotateLeft();
   // Wavefunctions
   for (std::size_t i = 0; i < Krylov.size(); ++i)
   {
      Krylov[i].RotateLeftExpand();
   }
   // matrix elements of < Krylov[i+1] | H | Krylov[i] >
   for (std::size_t i = 0; i < kn1_H_kn.size(); ++i)
   {
      SuperblockStackRotateLeft(kn1_H_kn[i], Krylov[i+1], H, Krylov[i]);
   }
   // matrix elements of <Krylov[i] | Krylov[j] >
   for (std::size_t i = 0; i < Krylov.size()-1; ++i)
   {
      for (std::size_t j = i+1; j < Krylov.size(); ++j)
      {
         TransformStackRotateLeft(ki_kj(i,j), Krylov[i], Krylov[j]);
      }
   }
   this->DebugCheckBasis();
}

void SplitKrylov::ExpandLeft()
{
   for (std::size_t i = 0; i < Krylov.size(); ++i)
   {
      Krylov[i].Left() = prod(Krylov[i].Left(), Krylov[i].Center());
      Krylov[i].Center() = ExpandBasis2(Krylov[i].Left());
   }
   // matrix elements of < Krylov[i+1] | H | Krylov[i] >
   for (std::size_t i = 0; i < kn1_H_kn.size(); ++i)
   {
      SuperblockStackUpdateLeft(kn1_H_kn[i], Krylov[i+1], H, Krylov[i]);
   }
   // matrix elements of <Krylov[i] | Krylov[j] >
   for (std::size_t i = 0; i < Krylov.size()-1; ++i)
   {
      for (std::size_t j = i+1; j < Krylov.size(); ++j)
      {
         TransformStackUpdateLeft(ki_kj(i,j), Krylov[i], Krylov[j]);
      }
   }
   this->DebugCheckBasis();
}

void SplitKrylov::ExpandRight()
{
   for (std::size_t i = 0; i < Krylov.size(); ++i)
   {
      Krylov[i].Right() = prod(Krylov[i].Center(), Krylov[i].Right());
      Krylov[i].Center() = ExpandBasis1(Krylov[i].Right());
   }
   // matrix elements of < Krylov[i+1] | H | Krylov[i] >
   for (std::size_t i = 0; i < kn1_H_kn.size(); ++i)
   {
      SuperblockStackUpdateRight(kn1_H_kn[i], Krylov[i+1], H, Krylov[i]);
   }
   // matrix elements of <Krylov[i] | Krylov[j] >
   for (std::size_t i = 0; i < Krylov.size()-1; ++i)
   {
      for (std::size_t j = i+1; j < Krylov.size(); ++j)
      {
         TransformStackUpdateRight(ki_kj(i,j), Krylov[i], Krylov[j]);
      }
   }
   this->DebugCheckBasis();
}

MatrixOperator
SplitKrylov::ConstuctWavefunctionFromRitz(Vector<std::complex<double> > const& x)
{
   // start with the Krylov[0] component
   MatrixOperator Result = x[0] * Krylov[0].Center();

   for (std::size_t i = 1; i < x.size(); +i)
   {
      Result += x[i] * triple_prod(ki_kj(0,i).Left(),
                                   Krylov[i].Center(),
                                   herm(ki_kj(0,i).Right()));
   }
   return Result;
}

MatrixOperator
SplitKrylov::ConstructLeftDensityMatrixFromRitz(Vector<std::complex<double> > const& x)
{
   using LinearAlgebra::norm_2_sq;
   // start with the (0,0) entry
   MatrixOperator Rho = norm_2_sq(x[0]) * scalar_prod(Krylov[0].Center(), herm(Krylov[0].Center()));
   // now the (0,j) and (j,0) entries
   for (std::size_t j = 1; j < x.size(); ++j)
   {
      MatrixOperator Next = x[0] * conj(x[j]) * scalar_prod(triple_prod(Krylov[0].Center(),
                                                                        ki_kj(0,j).Left(),
                                                                        herm(Krylov[j].Center())),
                                                            herm(ki_kj(0,j).Left()));
      Rho += Next + adjoint(Next);
   }
   // Now the entries (i,j) and (j,i), for j >= i > 1
   for (std::size_t i = 1; i < x.size(); ++i)
   {
      // the (i,i) entry
      Rho += x[i] * conj(x[i]) * triple_prod(ki_kj(0,i).Left(),
                                             scalar_prod(Krylov[i].Center(),
                                                         herm(Krylov[i].Center())),
                                             herm(ki_kj(0,i).Left()));
      // (i,j) for j > i
      for (std::size_t j = i+1; j < x.size(); ++j)
      {
         MatrixOperator Next = x[i] * conj(x[0]) * triple_prod(ki_kj(0,i).Left(),
                                                                triple_prod(Krylov[i].Center(),
                                                                            ki_kj(i,j).Left(),
                                                                            herm(Krylov[j].Center())),
                                                                herm(ki_kj(0,j).Left()));
         Rho += Next + adjoint(Next);
      }
   }

   return Rho;
}

MatrixOperator
SplitKrylov::ConstructRightDensityMatrixFromRitz(Vector<std::complex<double> > const& x)
{
   MatrixOperator Rho;
   // now the (0,j) and (j,0) entries
   for (std::size_t j = 1; j < x.size(); ++j)
   {
      Rho += x[0] * conj(x[j]) * scalar_prod(triple_prod(herm(Krylov[0].Center()),
                                                         ki_kj(0,j).Right(),
                                                         Krylov[j].Center()),
                                             herm(ki_kj(0,j).Right()));
   }
   // Now the entries (i,j) and (j,i), for j >= i > 1
   for (std::size_t i = 1; i < x.size(); ++i)
   {
      // (i,j) for j > i
      for (std::size_t j = i+1; j < x.size(); ++j)
      {
         Rho += x[i] * conj(x[0]) * triple_prod(herm(ki_kj(0,i).Right()),
                                                triple_prod(herm(Krylov[i].Center()),
                                                            ki_kj(i,j).Right(),
                                                            Krylov[j].Center()),
                                                ki_kj(0,j).Right());
      }
   }

   // The (j,i) components are now the hermitian conjugate
   Rho += adjoint(Rho);

   // Now (0,0) entry
   Rho += LinearAlgebra::norm_2_sq(x[0]) * scalar_prod(herm(Krylov[0].Center()), Krylov[0].Center());
   // the (i,i) entries
   for (std::size_t i = 1; i < x.size(); ++i)
   {
      Rho += x[i] * conj(x[i]) * triple_prod(herm(ki_kj(0,i).Right()),
                                             scalar_prod(herm(Krylov[i].Center()),
                                                         Krylov[i].Center()),
                                             ki_kj(0,i).Right());
   }

   return Rho;
}

void SplitKrylov::ConstructKrylovBasis()
{
   // We assume that Krylov[0] is normalized
   DEBUG_CHECK(equal(norm_frob(Krylov[0].Center()), 1.0, 1E-5))(norm_frob(Krylov[0].Center()));

   sub_H = Matrix<std::complex<double> >(Krylov.size()-1, Krylov.size()-1, 0.0);

   for (std::size_t i = 0; i < Krylov.size()-1; ++i)
   {
      Krylov[i+1].Center() = operator_prod(conj(H.Center()),
                                           kn1_H_kn[i].Left(),
                                           Krylov[i].Center(),
                                           herm(kn1_H_kn[i].Right()));
#if 0
      Vector<MatrixOperator> kj_in_basis_i(i+1);
      Matrix<std::complex<double> > P(i+1, i+1);
      for (std::size_t j = 0; j <= i; ++j)
      {
         kj_in_basis_i[j] = triple_prod(herm(ki_kj(j,i+1).Left()),
                                        Krylov[j].Center(),
                                        ki_kj(j,i+1).Right());

         for (int p = 0; p <= j; ++p)
         {
            P(p,j) =

         sub_H(i,j) = inner_prod(kj_in_this_basis, Krylov[i+1].Center());
         sub_H(j,i) = conj(sub_H(i,j));
         Krylov[i+1].Center() -= (sub_H(i,j) / norm_frob_sq(kj_in_this_basis)) * kj_in_this_basis;
         //TRACE(sub_H(i,j))(inner_prod(kj_in_this_basis, Krylov[i+1].Center()))(i)(j);
      }
#endif
      // normalize the Krylov vector
      Krylov[i+1] *= (1.0 / norm_frob(Krylov[i+1]));
   }
}

void
SplitKrylov::TruncateLeft(double MixFactor, MatrixOperator Rho_k0L,
                          MatrixOperator Rho_k0R, int MaxStates)
{
   this->TruncateCommon(MixFactor, Rho_k0L, Rho_k0R, MaxStates, -1);
}

void
SplitKrylov::TruncateRight(double MixFactor, MatrixOperator Rho_k0L,
                          MatrixOperator Rho_k0R, int MaxStates)
{
   this->TruncateCommon(MixFactor, Rho_k0L, Rho_k0R, MaxStates, 1);
}

void
SplitKrylov::TruncateCommon(double MixFactor, MatrixOperator Rho_k0L,
                            MatrixOperator Rho_k0R, int MaxStates, int Direction)
{
   // Firstly, we construct the density operators
   std::vector<MatrixOperator> RhoL(Krylov.size());
   std::vector<MatrixOperator> RhoR(Krylov.size());

   RhoL[0] = MixFactor * Rho_k0L;
   RhoR[0] = MixFactor * Rho_k0R;
   if (MixFactor < 1)
   {
      RhoL[0] += (1.0 - MixFactor) * scalar_prod(Krylov[0].Center(), herm(Krylov[0].Center()));
      RhoR[0] += (1.0 - MixFactor) * scalar_prod(herm(Krylov[0].Center()), Krylov[0].Center());
   }

   // the rest of the Krylov vectors
   for (std::size_t i = 1; i < Krylov.size(); ++i)
   {
      if (MixFactor > 0)
      {
         RhoL[i] = MixFactor *
            operator_prod(trace_prod(prod(kn1_H_kn[i-1].Right(), RhoR[i-1]),
                                     herm(kn1_H_kn[i-1].Right())),
                          kn1_H_kn[i-1].Left(),
                          RhoL[i-1],
                          herm(kn1_H_kn[i-1].Left()));
         RhoR[i] = MixFactor *
            operator_prod(trace_prod(prod(kn1_H_kn[i-1].Left(), RhoL[i-1]),
                                     herm(kn1_H_kn[i-1].Left())),
                          kn1_H_kn[i-1].Right(),
                          RhoR[i-1],
                          herm(kn1_H_kn[i-1].Right()));

         if (MixFactor < 1)
         {
            RhoL[i] += (1.0-MixFactor) * scalar_prod(Krylov[i].Center(), herm(Krylov[i].Center()));
            RhoR[i] += (1.0-MixFactor) * scalar_prod(herm(Krylov[i].Center()), Krylov[i].Center());
         }
      }
      else
      {
         // MixFactor == 0, no need for the other rho
         if (Direction == -1)
            RhoL[i] = scalar_prod(Krylov[i].Center(), herm(Krylov[i].Center()));
         else
            RhoR[i] = scalar_prod(herm(Krylov[i].Center()), Krylov[i].Center());
      }
   }

   // Now apply the truncations
   for (std::size_t i = 0; i < Krylov.size(); ++i)
   {
      if (Direction == 1)
      {
         MatrixOperator U = this->ConstructRightTruncator(i, RhoR[i], MaxStates);
         this->TruncateKrylovRight(i, U);
      }
      else
      {
         CHECK_EQUAL(Direction, -1);
         MatrixOperator U = this->ConstructLeftTruncator(i, RhoL[i], MaxStates);
         this->TruncateKrylovLeft(i, U);
      }
   }
   this->DebugCheckBasis();
}

MatrixOperator SplitKrylov::ConstructLeftTruncator(int i, MatrixOperator const& Rho, int MaxStates)
{
   DensityMatrix<MatrixOperator> DM(Rho);
   DensityMatrix<MatrixOperator>::const_iterator E = DM.begin();
   int Count = 0;
   while (E != DM.end() && Count < MaxStates)
   {
      ++E; ++Count;
   }
   return DM.ConstructTruncator(DM.begin(), E);
}

MatrixOperator SplitKrylov::ConstructRightTruncator(int i, MatrixOperator const& Rho, int MaxStates)
{
   DensityMatrix<MatrixOperator> DM(Rho);
   DensityMatrix<MatrixOperator>::const_iterator E = DM.begin();
   int Count = 0;
   while (E != DM.end() && Count < MaxStates)
   {
      ++E; ++Count;
   }
   return DM.ConstructTruncator(DM.begin(), E);
}

void SplitKrylov::TruncateKrylovLeft(std::size_t i, MatrixOperator const& U)
{
   this->DebugCheckBasis();
   Krylov[i].Center() = prod(U, Krylov[i].Center(), Ident);
   Krylov[i].Left() = prod(Krylov[i].Left(), herm(U));
   if (i < kn1_H_kn.size())
      kn1_H_kn[i].Left() = prod(kn1_H_kn[i].Left(), herm(U));
   if (i > 0)
      kn1_H_kn[i-1].Left() = prod(U, kn1_H_kn[i-1].Left());
   for (std::size_t k = 0; k < i; ++k)
   {
      ki_kj(k,i).Left() = ki_kj(k,i).Left(), herm(U);
   }
   for (std::size_t j = i+1; j < Krylov.size(); ++j)
   {
      ki_kj(i,j).Left() = prod(U, ki_kj(i,j).Left(), Ident);
   }
   this->DebugCheckBasis();
}

void SplitKrylov::TruncateKrylovRight(std::size_t i, MatrixOperator const& U)
{
   this->DebugCheckBasis();
   Krylov[i].Center() = Krylov[i].Center() * herm(U);
   Krylov[i].Right() = prod(U, Krylov[i].Right());
   if (i < kn1_H_kn.size())
      kn1_H_kn[i].Right() = prod(kn1_H_kn[i].Right(), herm(U));
   if (i > 0)
      kn1_H_kn[i-1].Right() = prod(U, kn1_H_kn[i-1].Right());
   for (std::size_t k = 0; k < i; ++k)
   {
      ki_kj(k,i).Right() = ki_kj(k,i).Right(), herm(U);
   }
   for (std::size_t j = i+1; j < Krylov.size(); ++j)
   {
      ki_kj(i,j).Right() = U * ki_kj(i,j).Right();
   }
   this->DebugCheckBasis();
}

void SplitKrylov::DebugCheckBasis() const
{
   for (std::size_t i = 0; i < Krylov.size(); ++i)
   {
      if (i != Krylov.size()-1)
      {
         DEBUG_CHECK_EQUAL(kn1_H_kn[i].Left().Basis2(), Krylov[i].Center().Basis1());
         DEBUG_CHECK_EQUAL(kn1_H_kn[i].Right().Basis2(), Krylov[i].Center().Basis2());
      }
      if (i != Krylov.size()-1)
      {
         DEBUG_CHECK_EQUAL(kn1_H_kn[i].Left().Basis1(), Krylov[i+1].Center().Basis1());
         DEBUG_CHECK_EQUAL(kn1_H_kn[i].Right().Basis1(), Krylov[i+1].Center().Basis2());
      }

      for (std::size_t j = i+1; j < Krylov.size(); ++j)
      {
         DEBUG_CHECK_EQUAL(ki_kj(i,j).Left().Basis2(), Krylov[j].Center().Basis1());
         DEBUG_CHECK_EQUAL(ki_kj(i,j).Left().Basis1(), Krylov[i].Center().Basis1());
         DEBUG_CHECK_EQUAL(ki_kj(i,j).Right().Basis2(), Krylov[j].Center().Basis2());
         DEBUG_CHECK_EQUAL(ki_kj(i,j).Right().Basis1(), Krylov[i].Center().Basis2());
      }
   }
}
