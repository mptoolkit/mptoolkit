// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/simplekrylov.cpp
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

#include "simplekrylov.h"
#include "linearalgebra/eigen.h"
#include "tensor/tensor_eigen.h"
#include "matrixproduct/wavefunc-utils.h"
#include "lanczos-exponential.h"
#include <boost/none.hpp>
#include <fstream>

using MessageLogger::Logger;
using MessageLogger::msg_log;

struct SuperblockMultiply
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator value_type;
   typedef MatrixOperator argument_type;

   SuperblockMultiply(SimpleOperator const& Ham_,
                      MPStateComponent const& Left_,
                      MPStateComponent const& Right_);

   MatrixOperator operator()(MatrixOperator const& Psi) const
   {
      return operator_prod(Ham, Left, Psi, herm(Right));
   }

   SimpleOperator Ham;
   MPStateComponent Left, Right;
};

inline
SuperblockMultiply::SuperblockMultiply(SimpleOperator const& Ham_,
                                       MPStateComponent const& Left_,
                                       MPStateComponent const& Right_)
   : Ham(Ham_), Left(Left_), Right(Right_)
{
}

SimpleKrylov::SimpleKrylov(CenterWavefunction const& Psi0_,
                           SplitOperator const& Op_,
                           std::complex<double> Timestep_,
                           CenterWavefunction const& Psi1_,
                           bool UsePsi1HPsi0_)
   : Psi0(Psi0_), Ham(Op_), Psi1(Psi1_), UsePsi1HPsi0(UsePsi1HPsi0_),
     Ident(Psi0_.GetSymmetryList()),
     Timestep(Timestep_),
     TotalSweepNumber(0)
{
   // construct the OpMatrix elements for the right hand side

   typedef MPStateComponent Component;

   // Shift H so the Center matrix lines up with that of Psi1
   while (Ham.LeftSize() > Psi1.LeftSize()) Ham.RotateLeft();
   while (Ham.RightSize() > Psi1.RightSize()) Ham.RotateRight();
   // Shift Psi0 so the Center matrix lines up with that of x
   while (Psi0.LeftSize() > Psi1.LeftSize()) Psi0.RotateLeft();
   while (Psi0.RightSize() > Psi1.RightSize()) Psi0.RotateRight();

   if (UsePsi1HPsi0)
      InitializeSuperblockStack(Psi1_H_Psi0, Psi1, Ham, Psi0);
   InitializeSuperblockStack(Psi1_H_Psi1, Psi1, Ham, Psi1);
   InitializeTransformStack(Psi1_Psi0, Psi1, Psi0);
}

SimpleKrylov::SimpleKrylov(CenterWavefunction const& Psi0_,
                           SplitOperator const& Op_,
                           std::complex<double> Timestep_,
                           SuperblockOperator const& Psi0_H_Psi0,
                           bool UsePsi1HPsi0_)
   : Psi0(Psi0_), Ham(Op_), Psi1(Psi0_), UsePsi1HPsi0(UsePsi1HPsi0_),
     Psi1_H_Psi0(Psi0_H_Psi0),
     Psi1_H_Psi1(Psi0_H_Psi0),
     Ident(Psi0_.GetSymmetryList()),
     Timestep(Timestep_)
{
   // Optimization: this should be a stack of identity operators
   InitializeTransformStack(Psi1_Psi0, Psi1, Psi0);
}

SimpleKrylov::~SimpleKrylov()
{
}

void
SimpleKrylov::SetupLogFiles(std::string const& Prefix, bool Truncate)
{
   std::ios_base::openmode Mode = Truncate ? (std::ios_base::out | std::ios_base::trunc)
      : (std::ios_base::out | std::ios_base::app);

   SweepLog = boost::shared_ptr<std::ofstream>(new std::ofstream((Prefix + ".sweep").c_str(), Mode));
   SweepLog->precision(12);
   StepLog = boost::shared_ptr<std::ofstream>(new std::ofstream((Prefix + ".step").c_str(), Mode));
   StepLog->precision(12);

   Logger("SweepLog").SetStream(*SweepLog);
   Logger("SweepLog").SetThreshold(10);
   Logger("StepLog").SetStream(*StepLog);
   Logger("StepLog").SetThreshold(10);

   msg_log(1, "StepLog") << "#1_LeftSize #2_RightSize #3_NumStates #4_Time #5_Variance #6_SolverVariance #7_Truncation #8_RealTruncation #9SolverTol #10_Entropy\n";
   msg_log(1, "SweepLog") << "#1_SweepNum #2_NumMultiplies #3_WallTime #4_CPUTime #5_Truncation #6_RealTruncation #7_Entropy #8_SolverVariance #9_Variance #10_StDevBound\n";
}

void SimpleKrylov::AdvanceTime()
{
   Psi0 = Psi1;
   // Psi1_H_Psi1 doesn't need modifying
   if (UsePsi1HPsi0)
      Psi1_H_Psi0 = Psi1_H_Psi1;
   InitializeTransformStack(Psi1_Psi0, Psi1, Psi0); // this will be identity operators
}

void SimpleKrylov::SetWavefunctionAttributes() const
{
   double RealTime = Psi0.Attributes()["Time"].get_or_default<double>(0.0);
   double Beta = Psi0.Attributes()["Beta"].get_or_default<double>(0.0);

   RealTime += Timestep.imag();
   Beta -= Timestep.real();

   if (RealTime != 0 || Psi0.Attributes().count("Time"))
      Psi1.AttributesMutable()["Time"] = RealTime;
   if (Beta != 0 || Psi0.Attributes().count("Beta"))
      Psi1.AttributesMutable()["Beta"] = Beta;
}

CenterWavefunction& SimpleKrylov::Wavefunction()
{
   this->SetWavefunctionAttributes();
   return Psi1;
}

CenterWavefunction const& SimpleKrylov::Wavefunction() const
{
   this->SetWavefunctionAttributes();
   return Psi1;
}

MatrixOperator const& SimpleKrylov::k0() const
{
   if (!k0Cache)
      k0Cache = triple_prod(Psi1_Psi0.Left(), Psi0.Center(), herm(Psi1_Psi0.Right()));
   return k0Cache.get();
}

void SimpleKrylov::Invalidatek0()
{
   k0Cache = boost::none;
}

double SimpleKrylov::Solve(int NumIterations, double ErrBound)
{
   OldPsi1Center = Psi1.Center();

   // Krylov vector k0 in the Psi1 basis
   MatrixOperator k0_ = this->k0();

   // this is what we would do if we used the matrix elements Psi1_H_Psi0
   // Krylov vector k1 in the Psi1 basis
   if (UsePsi1HPsi0)
   {
      MatrixOperator k1 = operator_prod(conj(Ham.Center()),
                                        Psi1_H_Psi0.Left(), Psi0.Center(), herm(Psi1_H_Psi0.Right()));
      Psi1.Center() = LanczosExponential(k0_, k1, SuperblockMultiply(conj(Ham.Center()),
                                                                     Psi1_H_Psi1.Left(),
                                                                     Psi1_H_Psi1.Right()),
                                         NumIterations,
                                         Timestep);
   }
   else
   {
      Psi1.Center() = LanczosExponential(k0_, SuperblockMultiply(conj(Ham.Center()),
                                                                 Psi1_H_Psi1.Left(),
                                                                 Psi1_H_Psi1.Right()),
                                         NumIterations,
                                         Timestep,
                                         ErrBound);
      IterationSolverVariance = norm_frob_sq(OldPsi1Center - Psi1.Center()) / (norm_frob(OldPsi1Center) * norm_frob(Psi1.Center()));
      IterationSolverTolerance = ErrBound;
   }

   Psi1.Center() *= 1.0 / norm_frob(Psi1.Center());

   IterationNumMultiplies += NumIterations;
   return NumIterations;
}

void SimpleKrylov::StartIteration()
{
   IterationNumMultiplies = 0;
   IterationNumStates = 0;
   IterationTruncation = 0;
   IterationRealTruncation = 0;
   IterationEntropy = 0;
   IterationSolverVariance = 0;
   IterationOverallVariance = 0;
   IterationSolverTolerance = 0;
   IterationStartTime = ProcControl::GetWallTime();
}

void SimpleKrylov::EndIteration()
{
   double IterationTime = ProcControl::GetWallTime() - IterationStartTime;

   msg_log(1, "StepLog") << this->LeftSize() << ' '
                         << this->RightSize() << ' '
                         << IterationNumStates << ' '
                         << IterationTime << ' '
                         << IterationOverallVariance << ' '
                         << IterationSolverVariance << ' '
                         << IterationTruncation << ' '
                         << IterationRealTruncation << ' '
                         << IterationEntropy << ' '
                         << '\n';

   ++SweepNumIterations;
   SweepNumMultiplies += IterationNumMultiplies;
   SweepSumStates += IterationNumStates;
   SweepMaxStates = std::max(SweepMaxStates, IterationNumStates);
   SweepTruncation += IterationTruncation;
   SweepEntropy = std::max(SweepEntropy, IterationEntropy);
   SweepSolverVariance += IterationSolverVariance;
   SweepOverallVariance += IterationOverallVariance;
   SweepStDevBound += std::sqrt(IterationOverallVariance);
   SweepRealTruncation += IterationRealTruncation;
}

void SimpleKrylov::StartSweep(bool)
{
   ++TotalSweepNumber;
   SweepNumIterations = 0;
   SweepSumStates = 0;
   SweepMaxStates = 0;
   SweepNumMultiplies = 0;
   SweepTruncation = 0;
   SweepRealTruncation = 0;
   SweepEntropy = 0;
   SweepSolverVariance = 0;
   SweepOverallVariance = 0;
   SweepStDevBound = 0;
   SweepStartTime = ProcControl::GetWallTime();
   SweepStartCPUTime = ProcControl::GetCPUTime();
}

void SimpleKrylov::EndSweep()
{
   double SweepTime = ProcControl::GetWallTime() - SweepStartTime;
   double SweepCPUTime = ProcControl::GetCPUTime() - SweepStartCPUTime;

   msg_log(1, "SweepLog") << TotalSweepNumber << ' '
                          << SweepNumMultiplies << ' '
                          << SweepTime << ' '
                          << SweepCPUTime << ' '
                          << SweepTruncation << ' '
                          << SweepRealTruncation << ' '
                          << SweepEntropy << ' '
                          << SweepSolverVariance << ' '
                          << SweepOverallVariance << ' '
                          << SweepStDevBound << ' '
                          << std::endl;

   msg_log(1, "StepLog") << std::flush;
}

void SimpleKrylov::ExpandLeft()
{
   // expand Psi1
   Psi1.Left() = prod(Psi1.Left(), Psi1.Center());
   Psi1.Center() = ExpandBasis2(Psi1.Left());

   // matrix elements
   if (UsePsi1HPsi0)
      SuperblockStackUpdateLeft(Psi1_H_Psi0, Psi1, Ham, Psi0);
   SuperblockStackUpdateLeft(Psi1_H_Psi1, Psi1, Ham, Psi1);
   TransformStackUpdateLeft(Psi1_Psi0, Psi1, Psi0);
   this->Invalidatek0();

   this->DebugCheckBasis();
}

void SimpleKrylov::ExpandRight()
{
   // expand Psi1
   Psi1.Right() = prod(Psi1.Center(), Psi1.Right());
   Psi1.Center() = ExpandBasis1(Psi1.Right());

   // matrix elements
   if (UsePsi1HPsi0)
      SuperblockStackUpdateRight(Psi1_H_Psi0, Psi1, Ham, Psi0);
   SuperblockStackUpdateRight(Psi1_H_Psi1, Psi1, Ham, Psi1);
   TransformStackUpdateRight(Psi1_Psi0, Psi1, Psi0);
   this->Invalidatek0();

   this->DebugCheckBasis();
}

void SimpleKrylov::ShiftLeft()
{
   Ham.RotateLeft();
   Psi0.RotateLeft();
   Psi1.RotateLeft();

   if (UsePsi1HPsi0)
      SuperblockStackRotateLeft(Psi1_H_Psi0, Psi1, Ham, Psi0);
   SuperblockStackRotateLeft(Psi1_H_Psi1, Psi1, Ham, Psi1);
   TransformStackRotateLeft(Psi1_Psi0, Psi1, Psi0);
   this->Invalidatek0();

   this->DebugCheckBasis();
}

void SimpleKrylov::ShiftRight()
{
   Ham.RotateRight();
   Psi0.RotateRight();
   Psi1.RotateRight();

   if (UsePsi1HPsi0)
      SuperblockStackRotateRight(Psi1_H_Psi0, Psi1, Ham, Psi0);
   SuperblockStackRotateRight(Psi1_H_Psi1, Psi1, Ham, Psi1);
   TransformStackRotateRight(Psi1_Psi0, Psi1, Psi0);
   this->Invalidatek0();

   this->DebugCheckBasis();
}

void SimpleKrylov::ShiftRightAndExpand()
{
   // Rotate the wavefunctions to the right, expand left basis of Psi1

   Ham.RotateRight();
   Psi0.RotateRight();
   Psi1.PushLeft(prod(Psi1.Center(), Psi1.Right()));
   Psi1.PopRight();
   Psi1.Center() = ExpandBasis2(Psi1.Left());

   if (UsePsi1HPsi0)
      SuperblockStackRotateRight(Psi1_H_Psi0, Psi1, Ham, Psi0);
   SuperblockStackRotateRight(Psi1_H_Psi1, Psi1, Ham, Psi1);
   TransformStackRotateRight(Psi1_Psi0, Psi1, Psi0);
   this->Invalidatek0();

   this->DebugCheckBasis();
}

void SimpleKrylov::ShiftLeftAndExpand()
{
   // Rotate the wavefunctions to the left, expand right basis of Psi1
   Ham.RotateLeft();
   Psi0.RotateLeft();
   Psi1.PushRight(prod(Psi1.Left(), Psi1.Center()));
   Psi1.PopLeft();
   Psi1.Center() = ExpandBasis1(Psi1.Right());

   if (UsePsi1HPsi0)
      SuperblockStackRotateLeft(Psi1_H_Psi0, Psi1, Ham, Psi0);
   SuperblockStackRotateLeft(Psi1_H_Psi1, Psi1, Ham, Psi1);
   TransformStackRotateLeft(Psi1_Psi0, Psi1, Psi0);
   this->Invalidatek0();

   this->DebugCheckBasis();
}

TruncationInfo SimpleKrylov::TruncateLeft(StatesInfo const& States, double MixFactor)
{
   MatrixOperator Rho1 = scalar_prod(Psi1.Center(), herm(Psi1.Center()));
   MatrixOperator Rho0 =  triple_prod(Psi1_Psi0.Left(),
                                      scalar_prod(Psi0.Center(), herm(Psi0.Center())),
                                      herm(Psi1_Psi0.Left()));
   //   MatrixOperator Rho0 = scalar_prod(this->k0(), herm(this->k0()));
   MatrixOperator Rho = 0.5 * (Rho1 + Rho0);

   if (MixFactor != 0)
   {
      MatrixOperator Mix;
      for (unsigned i = 0; i < Psi1_H_Psi1.Left().size(); ++i)
      {
         //         if (i != 0 && i !=  Psi1_H_Psi1.Left().size()-1)
         {
            double Normalizer = norm_frob_sq(Psi1_H_Psi1.Left()[i]) / Rho.Basis1().total_degree();
            MatrixOperator Correction = triple_prod(Psi1_H_Psi1.Left()[i],
                                                    Rho,
                                                    herm(Psi1_H_Psi1.Left()[i]));
            Correction *= 1.0 / Normalizer;
            Mix += Correction;
         }
      }
      Rho += MixFactor * Mix;
   }

   DensityMatrix<MatrixOperator> DM(Rho);
   TruncationInfo Info;
   MatrixOperator U = DM.ConstructTruncator(DM.begin(),
                                            TruncateFixTruncationErrorRelative(DM.begin(), DM.end(), States, Info));
   IterationNumStates = Info.KeptStates();
   IterationTruncation += Info.TruncationError();
   IterationEntropy = std::max(IterationEntropy, Info.KeptEntropy());

   // truncate the wavefunction
   Psi1.Center() = prod(U, Psi1.Center(), Ident);
   Psi1.Left() = prod(Psi1.Left(), herm(U));
   //Psi1.Normalize();

   // truncate the Ham matrix
   Psi1_H_Psi1.Left() = triple_prod(U, Psi1_H_Psi1.Left(), herm(U));
   if (UsePsi1HPsi0)
      Psi1_H_Psi0.Left() = prod(U, Psi1_H_Psi0.Left());

   // truncate the projector
   Psi1_Psi0.Left() = prod(U, Psi1_Psi0.Left(), Ident);

   // After the change in basis, our k0 is no longer valid
   this->Invalidatek0();

   // real truncation error
   IterationRealTruncation = std::max(0.0, 1.0 - norm_frob_sq(Psi1.Center()));

   // the truncation contributes to the wavefunction variance
   //   IterationVariance += 1.0 - norm_frob_sq(Psi1.Center());

   IterationOverallVariance = norm_frob_sq(OldPsi1Center - herm(U)*Psi1.Center());

   this->DebugCheckBasis();

  return Info;
}

TruncationInfo SimpleKrylov::TruncateRight(StatesInfo const& States, double MixFactor)
{
   MatrixOperator Rho1 = scalar_prod(herm(Psi1.Center()), Psi1.Center());
   MatrixOperator Rho0 =  triple_prod(Psi1_Psi0.Right(),
                                      scalar_prod(herm(Psi0.Center()), Psi0.Center()),
                                      herm(Psi1_Psi0.Right()));
      //scalar_prod(herm(this->k0()), this->k0());
   MatrixOperator Rho = 0.5 * (Rho1 + Rho0);

   if (MixFactor != 0)
   {
      MatrixOperator Mix;
      for (unsigned i = 0; i < Psi1_H_Psi1.Right().size(); ++i)
      {
         //         if (i != 0 && i !=  Psi1_H_Psi1.Right().size()-1)
         {
            double Normalizer = norm_frob_sq(Psi1_H_Psi1.Right()[i]) / Rho.Basis1().total_degree();
            MatrixOperator Correction = triple_prod(Psi1_H_Psi1.Right()[i],
                                                    Rho,
                                                    herm(Psi1_H_Psi1.Right()[i]));
            Correction *= 1.0 / Normalizer;
            Mix += Correction;
         }
      }

      Rho += MixFactor * Mix;
   }

   DensityMatrix<MatrixOperator> DM(Rho);
   TruncationInfo Info;
   MatrixOperator U = DM.ConstructTruncator(DM.begin(),
                                            TruncateFixTruncationErrorRelative(DM.begin(), DM.end(), States, Info));

   IterationNumStates = Info.KeptStates();
   IterationTruncation += Info.TruncationError();
   IterationEntropy = std::max(IterationEntropy, Info.KeptEntropy());

   // truncate the wavefunction
   Psi1.Center() = Psi1.Center() * herm(U);
   Psi1.Right() = prod(U, Psi1.Right());
   //Psi1.Normalize();

   // truncate the Ham matrix
   Psi1_H_Psi1.Right() = triple_prod(U, Psi1_H_Psi1.Right(), herm(U));
   if (UsePsi1HPsi0)
      Psi1_H_Psi0.Right() = prod(U, Psi1_H_Psi0.Right());

   // truncate the projector
   Psi1_Psi0.Right() = prod(U, Psi1_Psi0.Right(), Ident);

   // After the change in basis, our k0 is no longer valid
   this->Invalidatek0();

   // real truncation error
   IterationRealTruncation = 1.0 - norm_frob_sq(Psi1.Center());

   // the truncation contributes to the wavefunction variance
   //   IterationVariance += 1.0 - norm_frob_sq(Psi1.Center());

   IterationOverallVariance = norm_frob_sq(OldPsi1Center - Psi1.Center()*U);

   this->DebugCheckBasis();

   return Info;
}

void SimpleKrylov::DebugCheckBasis() const
{
   DEBUG_CHECK_EQUAL(Psi1_H_Psi1.Left().Basis1(), Psi1.Center().Basis1());
   DEBUG_CHECK_EQUAL(Psi1_H_Psi1.Right().Basis1(), Psi1.Center().Basis2());

   DEBUG_CHECK_EQUAL(Psi1_H_Psi1.Left().Basis2(), Psi1.Center().Basis1());
   DEBUG_CHECK_EQUAL(Psi1_H_Psi1.Right().Basis2(), Psi1.Center().Basis2());

   DEBUG_CHECK_EQUAL(Psi1_Psi0.Left().Basis1(), Psi1.Center().Basis1());
   DEBUG_CHECK_EQUAL(Psi1_Psi0.Right().Basis1(), Psi1.Center().Basis2());

   DEBUG_CHECK_EQUAL(Psi1_Psi0.Left().Basis2(), Psi0.Center().Basis1());
   DEBUG_CHECK_EQUAL(Psi1_Psi0.Right().Basis2(), Psi0.Center().Basis2());
}

PStream::opstream& operator<<(PStream::opstream& out, SimpleKrylov const& d)
{
   return out << d.Psi0
              << d.Ham
              << d.Psi1
              << d.Psi1_H_Psi0
              << d.Psi1_H_Psi1
              << d.Psi1_Psi0
              << d.Ident
              << d.Timestep

              << d.TotalSweepNumber

              << d.SweepNumIterations
              << d.SweepSumStates
              << d.SweepMaxStates
              << d.SweepNumMultiplies
              << d.SweepTruncation
              << d.SweepEntropy
              << d.SweepStartTime

              << d.IterationNumMultiplies
              << d.IterationNumStates
              << d.IterationTruncation
              << d.IterationEntropy;
}

PStream::ipstream& operator>>(PStream::ipstream& in, SimpleKrylov& d)
{
   return in >> d.Psi0
             >> d.Ham
             >> d.Psi1
             >> d.Psi1_H_Psi0
             >> d.Psi1_H_Psi1
             >> d.Psi1_Psi0
             >> d.Ident
             >> d.Timestep

             >> d.TotalSweepNumber

             >> d.SweepNumIterations
             >> d.SweepSumStates
             >> d.SweepMaxStates
             >> d.SweepNumMultiplies
             >> d.SweepTruncation
             >> d.SweepEntropy
             >> d.SweepStartTime

             >> d.IterationNumMultiplies
             >> d.IterationNumStates
             >> d.IterationTruncation
             >> d.IterationEntropy;
}
