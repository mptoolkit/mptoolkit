// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/dmrg2save.h
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

#if !defined(DMRG2_H_HHVERTHYEHYYOIUDLE89P)
#define DMRG2_H_HHVERTHYEHYYOIUDLE89P

#include "matrixproduct/mpstate.h"
#include "matrixproduct/mpoperator.h"
#include "matrixproduct/operatorstack.h"
#include "siteoperator/sitebasis.h"
#include "common/conflist.h"
#include <fstream>
#include <boost/shared_ptr.hpp>

typedef OperatorStack<MPMatrix> SuperblockOperator;

struct SuperblockMultiply
{
   typedef MPStateComponent result_type;
   typedef MPStateComponent argument_type;

   SuperblockMultiply(MPOpComponent const& H_, MPMatrix const& E_, MPMatrix const& F_)
      : H(H_), E(E_), F(F_) {}

   MPStateComponent operator()(MPStateComponent const& Psi) const
   {
      return operator_prod(H, E, Psi, herm(F));
   }

   MPOpComponent const& H;
   MPMatrix const& E;
   MPMatrix const& F;
};

struct SaveInfo
{
   MPStateComponent AFull;
   MPStateComponent HFull;
   DiagonalProjection UKeep;
   DiagonalProjection UEnv;
};

// When sweeping to the right, Psi is (KeptBasis, FullBasis)
// Env is (EnvBasis, FullBasis), RhoEnv is (EnvBasis, EnvBasis)
struct StateInfo
{
   MatrixOperator Psi;
   MatrixOperator Env;
   MatrixOperator RhoEnv;

   MPStateComponent HFull;
   DiagonalProjection UKeep;
   DiagonalProjection UEnv;
};

struct DMRG
{
   typedef OperatorStack<MPStateComponent> SuperblockOperator;
   typedef OperatorStack<MatrixOperator>   TransformOperator;

   DMRG() {}

   DMRG(MPWavefunction const& Psi_, MPOperator const& Ham_);

   void ShiftRight();
   void ShiftLeft();

   void TruncateLeft(int MinStates, int MaxStates, double MinTrunc, double CFactor);
   void TruncateRight(int MinStates, int MaxStates, double MinTrunc, double CFactor);

   void DoIterationRight();

   void debug_check_structure();

   MPOperator Ham;

   MPOperator::const_iterator Hi;

   OperatorStack<SaveInfo> Save;
   StateInfo State;

   int NumIterations;  // number of Lanczos iterations

   int MinStates;
   int MaxStates;
   double MinTrunc;
};

double Solve(MPStateComponent& Psi, MPStateComponent const& HLeft,
             MPOpComponent const& H, MPStateComponent const& HRight,
             int NumIter)
{
   std::vector<MPStateComponent> OrthoSet;
   double Energy = Lanczos(PsiCurrent, SuperblockMultiply(H, HLeft, HRight),
                           NumIter, OrthoSet);
   return Energy;
}

void DMRG::DoIterationRight()
{
   // Wavefunction and environment
   MPStateComponent Psi = State.Psi * Save.Right().AFull;
   MPStateComponent Env = State.Env * Save.Right().AFull;

   // Hamiltonian matrices
   MPStateComponent HLeft = triple_prod(State.UKeep, State.HFull, herm(State.UKeep));
   MPStateComponent HRight = Save.Right().HFull;

   // solve
   Solve(Psi, HLeft, *Hi, HRight, NumIterations);

   // Prepare for truncation
   MatrixOperator C = ExpandBasis2(Psi);

   // switch to the basis of DM eigenstates
   MatrixOperator Rho = scalar_prod(C, herm(C));
   DensityMatrix<MatrixOperator> DM(Rho);
   TruncationInfo Info;
   MatrixOperator U = DM.ConstructTruncator(DM.begin(), DM.end(), Info);
   C = prod(U, C);
   Psi = prod(Psi, herm(U));

   // Truncate into kept and discarded states
   DiagonalProjection UKeep, UDiscard;
   std::tie(UKeep, UDiscard) = basis_truncation(U.Basis2(), DM.begin(),
                                                  TruncateFixTruncationError(DM.begin(),
                                                                             DM.end(),
                                                                             MinStates,
                                                                             MaxStates,
                                                                             MinTrunc));

   // Discarded states
   MatrixOperator PsiDiscard = prod(UDiscard, C);

   // Kept states
   Psi = prod(Psi, herm(UKeep));
   C = prod(UKeep, C);

   // prepare for truncation of the environment
   MatrixOperator CEnv = ExpandBasis2(Env);

   // Merge the discarded basis with the expanded environment basis
   DiagonalProjection UComDiscard, UComEnv;
   std::tie(UComDiscard, UComEnv) = basis_sum(PsiDiscard.Basis1(), CEnv.Basis1());

   // Merge the discarded wavefunction with the environment, to get the remainder.
   MatrixOperator Rem = herm(UComDiscard)*PsiDiscard*UPsi + herm(UComEnv)*CEnv;

   // density matrix for the remainder
   MatrixOperator RhoRem = scalar_prod(Rem, herm(Rem));

   // Truncate the remainder, this becomes the new environment
   DensityMatrix<MatrixOperator> RemDM(RhoRem);
   TruncationInfo RemInfo;
   MatrixOperator RemTrunc = RemDM.ConstructTruncator(RemDM.begin(),
                                                      TruncateFixTruncationError(RemDM.begin(),
                                                                                 RemDM.end(),
                                                                                 MinStates,
                                                                                 MaxStates,
                                                                                 MinTrunc),
                                                      RemInfo);

   // Truncate the environment
   Env = prod(Env, herm(RemTrunc));
   CEnv = prod(RemTrunc, CEnv);

   // Merge the new environment with the kept states, to give the new full basis
   DiagonalProjection NewUKeep, NewUEnv;
   std::tie(NewUKeep, NewUEnv) = basis_sum(UKeep.Basis1(), RemTrunc.Basis1());

   // Get the full A matrix and make C and CEnv map from the full basis to the
   // kept and environment states
   MPStateComponent AFull = prod(Psi, NewUKeep) + prod(Env, NewUEnv);

   C = prod(herm(NewUKeep), C);
   CEnv = prod(herm(NewUEnv), CEnv);

   // Update Hamiltonian
   MPStateComponent HFull = operator_prod(*Hi, AFull, State.HFull, herm(AFull));

   // new SaveInfo
   SaveInfo Left;
   Left.AFull = AFull;
   Left.HFull = HFull;
   Left.UKeep = NewUKeep;
   Left.UEnv = NewUEnv;

   Save.PopRight();
   Save.PushLeft(Left);

   // new StateInfo
   State.Psi = C;
   State.Env = CEnv;
   State.HFull = HFull;
   State.UKeep = NewUKeep;
   State.UEnv = NewUEnv;

   // Shift the Hamiltonian
   ++Hi;
}

PStream::opstream& operator<<(PStream::opstream& out, DMRG const& d);
PStream::ipstream& operator>>(PStream::ipstream& out, DMRG& d);

DMRG::DMRG(MPWavefunction const& Psi_, MPOperator const& Ham_)
: Ham(Ham_)
{
   MPWavefunction::const_iterator Pi = Psi_.end();
   Hi = Ham.end();
   H.PushRight(make_vacuum_state(Psi_.GetSymmetryList()));
   while (Pi != Psi_.begin())
   {
      --Pi; --Hi;
      Psi.PushRight(*Pi, MatrixOperator::make_identity(Pi->Basis1()));
      H.PushRight(operator_product(*Hi, Psi.Right().A, E, herm(Psi.Right().A)));
   }
   CHECK_EQUAL(Hi, Ham.begin());
   H.PushLeft(make_vacuum_state(Psi_.GetSymmetryList()));
   PsiCurrent = Psi.Right().A;
   Psi.PopRight();
   H.PopRight();
}



double DMRG::Solve(int NumIter)
{
   return Energy;
}

void DMRG::ShiftRight()
{
   PsiMatrices.Right() = prod(Psi, PsiMatrices.Right());
   Psi = ExpandBasis1(PsiMatrices.Right());

   EnvMatrices.Right() = prod(Env, EnvMatrices.Right());
   OffDiagMatrices.Right() = prod(Env, OffDiagMatrices.Right());
   Env = ExpandBasis1(EnvMatrices.Right()) + ExpandBasis1(OffDiagMatrices.Right());

   // extend the right basis to the full basis
   DiagonalProjection UPsi, UEnv;
   std::tie(UPsi, UEnv) = basis_sum(Psi.Basis2(), Env.Basis2());

   PsiMatrices.Right() = herm(UPsi)*PsiMatrices.Right();
   Psi = prod(Psi, UPsi);

   EnvMatrices.Right() = herm(UEnv)*EnvMatrices.Right();
   OffDiagMatrices.Right() = herm(UEnv)*OffDiagMatrices.Right();
   Env = prod(Env, UEnv);

   // H operator
   HEnv.PopRight();
   HEnv.PushLeft(triple_prod(OffDiagMatrices.Left(), HPsi.Left(), herm(OffDiagMatrices.Left())));
   HPsi.PopRight();
   HPsi.PushLeft(triple_prod(PsiMatrices.Left(), HPsi.Left(), herm(PsiMatrices.Left())));
   HPsiEnv.PopRight();
   HPsiEnv.PushLeft(triple_prod(PsiMatrices.Left(), HPsiEnv.Left(), herm(EnvMatrices.Left())));
   HEnvPsi.PopRight();
   HEnv.Left() += triple_prod(EnvMatrices.Left(), HEnv.Left(), herm(EnvMatrices.Left()));
}

void DMRG::ShiftLeft()
{
   PsiMatrices.Left() = prod(PsiMatrices.Left(), Psi);
   Psi = ExpandBasis2(PsiMatrices.Left());

   EnvMatrices.Left() = prod(EnvMatrices.Right(), Env);
   OffDiagMatrices.Left() = prod(OffDiagMatrices.Left(), Env);
   Env = ExpandBasis2(EnvMatrices.Left()) + ExpandBasis2(OffDiagMatrices.Left());

   // Extend the left basis to the full basis
   DiagonalProjection UPsi, UEnv;
   std::tie(UPsi, UEnv) = basis_sum(Psi.Basis2(), Env.Basis2());

   PsiMatrices.Right() = herm(UPsi)*PsiMatrices.Right();
   Psi = prod(Psi, UPsi);

   EnvMatrices.Right() = herm(UEnv)*EnvMatrices.Right();
   OffDiagMatrices.Right() = herm(UEnv)*OffDiagMatrices.Right();
   Env = prod(Env, UEnv);

   // H operator
   HEnv.PopRight();
   HEnv.Left() = HEnv.Left() + triple_prod(OffDiagMatrices.Left(), HPsi.Left(), herm(OffDiagMatrices.Left()));
   HPsi.PopRight();
   HPsi.PushLeft(triple_prod(PsiMatrices.Left(), HPsi.Left(), herm(PsiMatrices.Left())));
   HEnv.PopRight();
   HEnv.PushLeft(triple_prod(EnvMatrices.Left(), HEnv.Left(), herm(EnvMatrices.Left())));
}

void DMRG::TruncateLeft(int MinStates, int MaxStates, double MinTrunc, double CFactor)
{
#if 0
   // We might need the single site density matrix later
   SimpleOperator SingleSiteDM;
   if (CFactor > 0)
   {
      MPStateComponent SS = prod(PsiMatrices.top().A, Psi);
      SingleSiteDM = trace_prod(SS, herm(SS));
   }
#endif

   MatrixOperator Rho = scalar_prod(Psi, herm(Psi));

   // switch to the basis of DM eigenstates
   DensityMatrix<MatrixOperator> DM(Rho);
   TruncationInfo Info;
   MatrixOperator U = DM.ConstructTruncator(DM.begin(), DM.end(), Info);
   Psi = prod(U, Psi);
   PsiMatrices.Left() = prod(PsiMatrices.Left(), herm(U));

   // Truncate into kept and discarded states
   DiagonalProjection UKeep, UDiscard;
   std::tie(UKeep, UDiscard) = basis_truncation(U.Basis2(), DM.begin(),
                                                  TruncateFixTruncationError(DM.begin(),
                                                                             DM.end(),
                                                                             MinStates,
                                                                             MaxStates,
                                                                             MinTrunc));

   // Discarded states
   MatrixOperator PsiDiscard = prod(UDiscard, Psi);

   // the off-diagonal A-matrix from the kept to discarded states

   // Merge the discarded basis with the expanded environment basis
   DiagonalProjection UComDiscard, UComEnv;
   std::tie(UComDiscard, UComEnv) = basis_sum(PsiDiscard.Basis1(), Env.Basis1());

   // Merge the discarded wavefunction with the environment, to get the remainder.
   // Transform the right basis to the full basis at the same time
   MatrixOperator Rem = herm(UComDiscard)*PsiDiscard*UPsi + herm(UComEnv)*Env*UEnv;

   // density matrix for the remainder
   MatrixOperator RhoRem = scalar_prod(Rem, herm(Rem));

#if 0
   if (CFactor > 0)
   {
      double RemWeight = trace(RhoRem);
      // environment density matrix in the new basis
      // Add some diagonal component, just in case the single site DM is singular
      SingleSiteDM += std::numeric_limits<double>::epsilon()
         * SimpleOperator::make_identity(SingleSiteDM.Basis1());
      SingleSiteDM *= 1.0 / trace(SingleSiteDM);
      RhoRem *= (1.0 - CFactor);
      RhoRem += CFactor * RemWeight * triple_prod(SingleSiteDM, herm(PsiEnv), RhoRem, PsiEnv);
   }
#endif

   // Truncate the remainder
   DensityMatrix<MatrixOperator> RemDM(RhoRem);
   TruncationInfo RemInfo;
   MatrixOperator RemTrunc = RemDM.ConstructTruncator(RemDM.begin(),
                                                      TruncateFixTruncationError(RemDM.begin(),
                                                                                 RemDM.end(),
                                                                                 MinStates,
                                                                                 MaxStates,
                                                                                 MinTrunc),
                                                      RemInfo);

   // The 'discarded' component of the remainder, forms the off-diagonal part
   MatrixOperator UDiscardEnv = RemTrunc * herm(UComDiscard);
   OffDiagMatrices.PushLeft(prod(Psi.Left(), herm(UDiscardEnv)));

   // The diagonal component of the remainder
   MatrixOperator UEnvEnv = RemTrunc * herm(UComEnv);
   EnvMatrices.Left() = prod(EnvMatrices.Left(), herm(UEnvEnv));

   // New environment; RemTrunc = UDiscardEnv \oplus UEnvEnv
   Env = prod(RemTrunc, Rho);

   // Truncate to Kept states
   Psi = prod(UKeep, Psi);
   PsiMatrices.Left() = prod(PsiMatrices.left(), herm(UKeep));

   this->debug_check_structure();
}

void debug_check_structure() const
{
   DEBUG_CHECK_EQUAL(PsiMatrices.Left().Basis2(), Psi.Basis1());
   DEBUG_CHECK_EQUAL(OffDiagMatrices.Left().Basis2(), Env.Basis1());
   DEBUG_CHECK_EQUAL(EnvMatrices.Left().Basis2(), Env.Basis1());
}


#endif
