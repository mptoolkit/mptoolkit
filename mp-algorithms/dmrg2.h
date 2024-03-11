// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp-algorithms/dmrg2.h
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

#if !defined(DMRG2_H_HHVERTHYEHYYOIUDLE89P)
#define DMRG2_H_HHVERTHYEHYYOIUDLE89P

#include "matrixproduct/mpoperator.h"
#include "matrixproduct/linearwavefunction.h"
#include "matrixproduct/leftrightstack.h"
#include "matrixproduct/vectorbasissum.h"
#include "siteoperator/sitebasis.h"
#include "lanczos-ortho.h"
#include "common/conflist.h"
#include <fstream>
#include <boost/shared_ptr.hpp>

//typedef OperatorStack<MPMatrix> SuperblockOperator;

struct SuperblockMultiply
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator argument_type;

   SuperblockMultiply(MPMatrix const& E_, MPMatrix const& F_)
      : E(E_), F(F_) {}

   MatrixOperator operator()(MatrixOperator const& C) const
   {
      return operator_prod(E, C, herm(F));
   }

   MPMatrix const& E;
   MPMatrix const& F;
};

// When sweeping to the right, HFull is the Hamiltonian for AFull.Basis1(),
// UKeep, UEnv are the the kept and environment projectors for AFull.Basis2().
//
// When sweeping to the left, HFull is the Hamiltonian for AFull.Basis2(),
// UKeep, UEnv are the kept and environment projectors for AFull.Basis1().
struct SaveInfo
{
   MPStateComponent AFull;
   MPStateComponent HFull;
   DiagonalProjection UKeep;
   DiagonalProjection UEnv;
};

PStream::opstream& operator<<(PStream::opstream& out, SaveInfo const& x);
PStream::ipstream& operator>>(PStream::ipstream& in, SaveInfo& x);

PStream::opstream& operator<<(PStream::opstream& out, SaveInfo const& x)
{
   return out << x.AFull << x.HFull << x.UKeep << x.UEnv;
}

PStream::ipstream& operator>>(PStream::ipstream& in, SaveInfo& x)
{
   return in >> x.AFull >> x.HFull >> x.UKeep >> x.UEnv;
}

// When sweeping to the right, Psi is (KeptBasis, FullBasis)
// Env is (EnvBasis, FullBasis), RhoEnv is (EnvBasis, EnvBasis).
// UKeep,UEnv are the projectors onto the kept and environment states
// for the Basis1.
//
// When sweeping to the left, Psi is (FullBasis, KeptBasis)
// Env is (FullBasis, EnvBasis), RhoEnv is (EnvBasis, EnvBasis).
struct StateInfo
{
   MatrixOperator Psi;
   MatrixOperator Env;
   MatrixOperator RhoEnv;

   MPStateComponent HFull;
   DiagonalProjection UKeep;
   DiagonalProjection UEnv;
   KeepListType KeepList;
};

struct DMRG
{
   //   typedef OperatorStack<MPStateComponent> SuperblockOperator;
   //   typedef OperatorStack<MatrixOperator>   TransformOperator;

   DMRG() {}

   DMRG(LinearWavefunction const& Psi_, MPOperator const& Ham_);

   void DoIterationMoveRight();
   void DoIterationMoveLeft();

   // Switches the State from moving right to moving left
   void RightEndpoint();

   // Switches the State from moving left to moving right
   void LeftEndpoint();

   double Solve(MPMatrix const& HLeft, MatrixOperator& C, MPMatrix const& HRight);

   // Constructs the projectors into the kept and discarded states,
   // for the given density operator.  Also uses State.KeepList.
   // A nice improvement here would be to promote KeepList states from the
   // environment, if possible.  But this is messy.
   std::pair<MatrixOperator, MatrixOperator>
   ConstructSplitTruncator(MatrixOperator const& Rho,
                           std::set<QuantumNumbers::QuantumNumber> const& AddedQN,
                           TruncationInfo& Info);

   void debug_check_structure() const;

   MPOperator Ham;

   MPOperator::const_iterator Hi;

   LeftRightStack<SaveInfo> Save;
   StateInfo State;

   int NumIterations;  // number of Lanczos iterations

   double MixFactor;
   double EnvMixFactor;
   double DiscardMixFactor;

   int MinStates;
   int MaxStates;
   double MinTrunc;

   int EnvMinStates;
   int EnvMaxStates;
   double EnvKeptWeight;
};

PStream::opstream& operator<<(PStream::opstream& out, DMRG const& d);
PStream::ipstream& operator>>(PStream::ipstream& out, DMRG& d);

double DMRG::Solve(MPMatrix const& HLeft, MatrixOperator& C, MPMatrix const& HRight)
{
   int NumIter = NumIterations;
   std::vector<MatrixOperator> OrthoSet;
   MatrixOperator CSave = C;
   double Energy = Lanczos(C, SuperblockMultiply(HLeft, HRight),
                           NumIter, OrthoSet);

   // sign of the eigenvector
   if (inner_prod(CSave, C).real() < 0)
      C *= -1.0;

   return Energy;
}

std::pair<MatrixOperator, MatrixOperator>
DMRG::ConstructSplitTruncator(MatrixOperator const& Rho,
                              std::set<QuantumNumbers::QuantumNumber> const& AddedQN,
                              TruncationInfo& Info)
{
   DensityMatrix<MatrixOperator> DM(Rho);
   DensityMatrix<MatrixOperator>::const_iterator DMPivot =
      TruncateFixTruncationErrorAbsolute(DM.begin(),
                                         DM.end(),
                                         MinStates,
                                         MaxStates,
                                         MinTrunc,
                                         Info);
   std::list<EigenInfo> KeepStates(DM.begin(), DMPivot);
   std::list<EigenInfo> DiscardStates(DMPivot, DM.end());
   UpdateKeepList(State.KeepList,
                  AddedQN,
                  DM.Basis(),
                  KeepStates,
                  DiscardStates,
                  Info);

   MatrixOperator UKeep = DM.ConstructTruncator(KeepStates.begin(), KeepStates.end());
   MatrixOperator UDiscard = DM.ConstructTruncator(DiscardStates.begin(), DiscardStates.end());
   return std::pair<MatrixOperator, MatrixOperator>(UKeep, UDiscard);
}

void DMRG::DoIterationMoveRight()
{
   // Wavefunction and environment
   MPStateComponent Psi = prod(State.Psi, Save.Right().AFull);
   MPStateComponent Env = prod(State.Env, Save.Right().AFull);

   // Hamiltonian matrices
   MPStateComponent HLeft = triple_prod(State.UKeep, State.HFull, herm(State.UKeep));
   MPStateComponent HRight = Save.Right().HFull;

   // It would be better to solve using Psi directly, but we don't have the
   // operator_prod overload yet
   MatrixOperator C = ExpandBasis2(Psi);
   HLeft = operator_prod(herm(*Hi), herm(Psi), HLeft, Psi);
   this->Solve(HLeft, C, HRight);

   // Density matrix
   MatrixOperator Rho = scalar_prod(C, herm(C));
   // White's density matrix mixing
   if (MixFactor > 0)
   {
      MatrixOperator RhoMix = operator_prod(HLeft, Rho, herm(HLeft));
      RhoMix *= MixFactor / trace(RhoMix);
      Rho *= 1.0 - MixFactor;
      Rho += RhoMix;
   }
   MatrixOperator UKeep, UDiscard;
   TruncationInfo Info;
   std::tie(UKeep, UDiscard) =
      this->ConstructSplitTruncator(Rho,
                                    QuantumNumbersInBasis(adjoint(Psi.SiteBasis())),
                                    Info);

   TRACE(Info.TruncationError());

   // Wavefunction for the discarded states
   MatrixOperator PsiDiscard = UDiscard*C;

   // prepare for truncation of the environment
   MatrixOperator CEnv = ExpandBasis2(Env);

   // Merge the discarded basis with the expanded environment basis
   DiagonalProjection UComDiscard, UComEnv;
   std::tie(UComDiscard, UComEnv) = basis_sum(PsiDiscard.Basis1(), CEnv.Basis1());

   // Merge the discarded wavefunction with the environment, to get the remainder.
   MatrixOperator Rem = herm(UComDiscard)*PsiDiscard + herm(UComEnv)*CEnv;

   // density matrix for the remainder
   MatrixOperator RhoRem = scalar_prod(Rem, herm(Rem));
   TRACE(trace(RhoRem));

   // Mix in the discarded states
   if (DiscardMixFactor > 0)
   {
      //MatrixOperator DiscardIdent = MatrixOperator::make_identity(UComDiscard.Basis1());
      //DiscardIdent *= DiscardMixFactor / trace(DiscardIdent);
      //DiscardIdent = triple_prod(herm(UComDiscard), DiscardIdent, UComDiscard);
      //RhoRem += DiscardIdent;
      MatrixOperator RhoDiscard = triple_prod(UDiscard, Rho, herm(UDiscard));
      RhoRem += triple_prod(herm(UComDiscard), RhoDiscard, UComDiscard);
   }
   if (EnvMixFactor > 0 && !State.RhoEnv.is_null())
   {
      // Mix in part of the direct product of the old environment density operator and the
      // single-site density operator
      MPStateComponent A = prod(Psi, C);
      SimpleOperator SingleSiteDM = trace_prod(A, herm(A));
      SingleSiteDM += 1E-10 * SimpleOperator::make_identity(SingleSiteDM.Basis1());
      MatrixOperator RhoEnvSite = operator_prod(SingleSiteDM, herm(Env), State.RhoEnv, Env);
      RhoEnvSite = triple_prod(herm(UComEnv), RhoEnvSite, UComEnv);
      RhoEnvSite *= EnvMixFactor * trace(RhoRem) / (trace(RhoEnvSite) + std::numeric_limits<double>::epsilon());
      RhoRem *= 1.0 - EnvMixFactor;
      RhoRem += RhoEnvSite;
   }

   // Truncate the remainder, this becomes the new environment
   DensityMatrix<MatrixOperator> RemDM(RhoRem);
   TruncationInfo RemInfo;
   MatrixOperator RemTrunc = RemDM.ConstructTruncator(RemDM.begin(),
             TruncateFixKeptWeight(RemDM.begin(),
                                   RemDM.end(),
                                   EnvMinStates,
                                   EnvMaxStates,
                                   EnvKeptWeight,
                                   RemInfo));

   // Truncate the environment
   Rem = RemTrunc*Rem;
   if (EnvMixFactor > 0)
      RhoRem = triple_prod(RemTrunc, RhoRem, herm(RemTrunc));

   // Split the remainder into the discarded and environment parts
   MatrixOperator RemDiscard = RemTrunc*herm(UComDiscard);
   MatrixOperator RemEnv = RemTrunc*herm(UComEnv);

   // Split Psi into the kept and discarded parts
   MPStateComponent PsiRem = prod(Psi, herm(UDiscard));
   PsiRem = prod(PsiRem, herm(RemDiscard));
   PsiRem = prod(herm(State.UKeep), PsiRem);
   // Kept states for Psi
   Psi = prod(Psi, herm(UKeep));
   C = UKeep*C;

   // Expand Env into the remainder basis, and the left full basis
   Env = prod(Env, herm(RemEnv));
   Env = prod(herm(State.UEnv), Env) + PsiRem;
   CEnv = Rem;

   // Merge the new environment with the kept states, to give the new full basis
   DiagonalProjection NewUKeep, NewUEnv;
   std::tie(NewUKeep, NewUEnv) = basis_sum(UKeep.Basis1(), RemTrunc.Basis1());

   // Get the full A matrix and make C and CEnv map from the full basis to the
   // kept and environment states
   MPStateComponent AFull = prod(herm(State.UKeep), prod(Psi, NewUKeep)) + prod(Env, NewUEnv);

   // new SaveInfo
   SaveInfo Left;
   Left.AFull = AFull;
   Left.HFull = State.HFull;
   Left.UKeep = NewUKeep;
   Left.UEnv = NewUEnv;

   Save.PushLeft(Left);
   Save.PopRight();

   // new StateInfo
   State.Psi = C;
   State.Env = CEnv;
   State.HFull = operator_prod(herm(*Hi), herm(AFull), State.HFull, AFull);
   State.UKeep = NewUKeep;
   State.UEnv = NewUEnv;
   if (EnvMixFactor > 0)
      State.RhoEnv = RhoRem;

   // Shift the Hamiltonian
   ++Hi;
}

void DMRG::DoIterationMoveLeft()
{
   // Shift the Hamiltonian
   --Hi;

   // Wavefunction and environment
   MPStateComponent Psi = prod(Save.Left().AFull, State.Psi);
   MPStateComponent Env = prod(Save.Left().AFull, State.Env);

   // Hamiltonian matrices
   MPStateComponent HRight = triple_prod(State.UKeep, State.HFull, herm(State.UKeep));
   MPStateComponent HLeft = Save.Left().HFull;

   // It would be better to solve using Psi directly, but we don't have the
   // operator_prod overload yet
   MatrixOperator C = ExpandBasis1(Psi);
   HRight = operator_prod(*Hi, Psi, HRight, herm(Psi));
   this->Solve(HLeft, C, HRight);

   // Density matrix
   MatrixOperator Rho = scalar_prod(herm(C), C);
   // White's density matrix mixing
   if (MixFactor > 0)
   {
      MatrixOperator RhoMix = operator_prod(HRight, Rho, herm(HRight));
      RhoMix *= MixFactor / trace(RhoMix);
      Rho *= 1.0 - MixFactor;
      Rho += RhoMix;
   }
   MatrixOperator UKeep, UDiscard;
   TruncationInfo Info;
   std::tie(UKeep, UDiscard) =
      this->ConstructSplitTruncator(Rho,
                                    QuantumNumbersInBasis(adjoint(Psi.SiteBasis())),
                                    Info);

   TRACE(Info.TruncationError());

   // Wavefunction for the discarded states
   MatrixOperator PsiDiscard = C*herm(UDiscard);

   // prepare for truncation of the environment
   MatrixOperator CEnv = ExpandBasis1(Env);

   // Merge the discarded basis with the expanded environment basis
   DiagonalProjection UComDiscard, UComEnv;
   std::tie(UComDiscard, UComEnv) = basis_sum(PsiDiscard.Basis2(), CEnv.Basis2());

   // Merge the discarded wavefunction with the environment, to get the remainder.
   MatrixOperator Rem = PsiDiscard*UComDiscard + CEnv*UComEnv;

   // density matrix for the remainder
   MatrixOperator RhoRem = scalar_prod(herm(Rem), Rem);
   TRACE(trace(RhoRem));

   // Mix in the discarded states
   if (DiscardMixFactor > 0)
   {
      //      MatrixOperator DiscardIdent = MatrixOperator::make_identity(UComDiscard.Basis1());
      //      DiscardIdent *= DiscardMixFactor / trace(DiscardIdent);
      //      DiscardIdent = triple_prod(herm(UComDiscard), DiscardIdent, UComDiscard);
      //      RhoRem += DiscardIdent;
      MatrixOperator RhoDiscard = triple_prod(UDiscard, Rho, herm(UDiscard));
      RhoRem += triple_prod(herm(UComDiscard), RhoDiscard, UComDiscard);
   }
   if (EnvMixFactor > 0 && !State.RhoEnv.is_null())
   {
      // Mix in part of the direct product of the old environment density operator and the
      // single-site density operator
      MPStateComponent A = prod(C, Psi);
      SimpleOperator SingleSiteDM = trace_prod(A, herm(A));
      SingleSiteDM += 1E-10 * SimpleOperator::make_identity(SingleSiteDM.Basis1());
      MatrixOperator RhoEnvSite = operator_prod(SingleSiteDM, Env, State.RhoEnv, herm(Env));
      RhoEnvSite = triple_prod(herm(UComEnv), RhoEnvSite, UComEnv);
      RhoEnvSite *= EnvMixFactor * trace(RhoRem) / (trace(RhoEnvSite) + std::numeric_limits<double>::epsilon());
      RhoRem *= 1.0 - EnvMixFactor;
      RhoRem += RhoEnvSite;
   }

   // Truncate the remainder, this becomes the new environment
   DensityMatrix<MatrixOperator> RemDM(RhoRem);
   TruncationInfo RemInfo;
   MatrixOperator RemTrunc = RemDM.ConstructTruncator(RemDM.begin(),
             TruncateFixKeptWeight(RemDM.begin(),
                                   RemDM.end(),
                                   EnvMinStates,
                                   EnvMaxStates,
                                   EnvKeptWeight,
                                   RemInfo));

   // Truncate the environment
   Rem = Rem*herm(RemTrunc);
   if (EnvMixFactor > 0)
      RhoRem = triple_prod(RemTrunc, RhoRem, herm(RemTrunc));

   // Split the remainder into the discarded and environment parts
   MatrixOperator RemDiscard = RemTrunc*herm(UComDiscard);
   MatrixOperator RemEnv = RemTrunc*herm(UComEnv);

   // Split Psi into the kept and discarded parts
   MPStateComponent PsiRem = prod(UDiscard, Psi);
   PsiRem = prod(RemDiscard, PsiRem);
   PsiRem = prod(PsiRem, State.UKeep);
   // Kept states for Psi
   Psi = prod(UKeep, Psi);
   C = C*herm(UKeep);

   // Expand Env into the remainder basis, and the left full basis
   Env = prod(RemEnv, Env);
   Env = prod(Env, State.UEnv) + PsiRem;
   CEnv = Rem;

   // Merge the new environment with the kept states, to give the new full basis
   DiagonalProjection NewUKeep, NewUEnv;
   std::tie(NewUKeep, NewUEnv) = basis_sum(UKeep.Basis1(), RemTrunc.Basis1());

   // Get the full A matrix and make C and CEnv map from the full basis to the
   // kept and environment states
   MPStateComponent AFull = prod(prod(herm(NewUKeep), Psi), State.UKeep) + prod(herm(NewUEnv), Env);

   // new SaveInfo
   SaveInfo Right;
   Right.AFull = AFull;
   Right.HFull = State.HFull;
   Right.UKeep = NewUKeep;
   Right.UEnv = NewUEnv;

   Save.PushRight(Right);
   Save.PopLeft();

   // new StateInfo
   State.Psi = C;
   State.Env = CEnv;
   State.HFull = operator_prod(*Hi, AFull, State.HFull, herm(AFull));
   State.UKeep = NewUKeep;
   State.UEnv = NewUEnv;
   if (EnvMixFactor > 0)
      State.RhoEnv = RhoRem;
}

void DMRG::RightEndpoint()
{
   // Expand the left basis for Psi to the full basis
   State.Psi = herm(State.UKeep) * State.Psi;
   // Normalize Psi
   double Norm = norm_frob(State.Psi);
   State.Psi *= 1.0 / Norm;
   TRACE(Norm);
   // Reset the Environment to zeron
   VectorBasis ZeroBasis(State.Psi.GetSymmetryList());
   State.Env = MatrixOperator(State.UEnv.Basis2(), ZeroBasis);
   State.UKeep = MatrixOperator::make_identity(State.Psi.Basis2());
   State.UEnv = MatrixOperator(ZeroBasis, State.Psi.Basis2());
   State.HFull = make_vacuum_state(State.Psi.GetSymmetryList());
   State.KeepList.clear();
   State.RhoEnv = MatrixOperator();
}

void DMRG::LeftEndpoint()
{
   // Expand the right basis for Psi to the full basis
   State.Psi = State.Psi * State.UKeep;
   // Normalize Psi
   double Norm = norm_frob(State.Psi);
   State.Psi *= 1.0 / Norm;
   TRACE(Norm);
   // Reset the Environment to zeron
   VectorBasis ZeroBasis(State.Psi.GetSymmetryList());
   State.Env = MatrixOperator(ZeroBasis, State.UEnv.Basis2());
   State.UKeep = MatrixOperator::make_identity(State.Psi.Basis1());
   State.UEnv = MatrixOperator(ZeroBasis, State.Psi.Basis1());
   State.HFull = make_vacuum_state(State.Psi.GetSymmetryList());
   State.KeepList.clear();
   State.RhoEnv = MatrixOperator();
}

DMRG::DMRG(LinearWavefunction const& Psi_, MPOperator const& Ham_)
: Ham(Ham_)
{
   LinearWavefunction::const_iterator Pi = Psi_.end();
   Hi = Ham.end();
   SaveInfo S;
   S.HFull = make_vacuum_state(Psi_.GetSymmetryList());
   S.UKeep = MatrixOperator::make_identity(S.HFull.Basis1()); // trivial
   VectorBasis ZeroBasis(S.HFull.GetSymmetryList());
   S.UEnv = MatrixOperator(ZeroBasis, S.HFull.Basis1()); // zero
   while (Pi != Psi_.begin())
   {
      --Pi; --Hi;
      S.AFull = *Pi;
      S.UKeep = MatrixOperator::make_identity(S.AFull.Basis1());
      S.UEnv = MatrixOperator(ZeroBasis, S.AFull.Basis1()); // zero
      Save.PushRight(S);
      S.HFull = operator_prod(*Hi, S.AFull, S.HFull, herm(S.AFull));
   }
   CHECK_EQUAL(Hi, Ham.begin());

   State.Psi = S.UKeep;
   State.Env = S.UEnv; // MatrixOperator(S.HFull.Basis1(), ZeroBasis);
   State.HFull = make_vacuum_state(Psi_.GetSymmetryList());
   State.UKeep = S.UKeep;
   State.UEnv = S.UEnv;
   std::set<QuantumNumbers::QuantumNumber> qs(QuantumNumbersInBasis(State.Psi.Basis1()));
   State.KeepList.insert(qs.begin(), qs.end());
}

void DMRG::debug_check_structure() const
{
   //   DEBUG_CHECK_EQUAL(PsiMatrices.Left().Basis2(), Psi.Basis1());
   //   DEBUG_CHECK_EQUAL(OffDiagMatrices.Left().Basis2(), Env.Basis1());
   //   DEBUG_CHECK_EQUAL(EnvMatrices.Left().Basis2(), Env.Basis1());
}


#endif
