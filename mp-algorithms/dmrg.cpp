// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/dmrg.cpp
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
// $Id$

#include "dmrg.h"
#include "lanczos-ortho.h"
#include "davidson.h"
#include "arnoldi.h"
#include "matrixproduct/density.h"
#include "pstream/optional.h"
#include <boost/optional.hpp>
#include <boost/none.hpp>
#include "linearalgebra/matrix_utility.h"
#include "matrixproduct/wavefunc-utils.h"
#include "common/statistics.h"
#include "common/environment.h"

using MessageLogger::Logger;
using MessageLogger::msg_log;

struct SuperblockMultiply
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator value_type;
   typedef MatrixOperator argument_type;

   SuperblockMultiply(MPStateComponent const& Left_,
                      MPStateComponent const& Right_);

   MatrixOperator operator()(MatrixOperator const& Psi) const
   {
      return operator_prod(Left, Psi, herm(Right));
   }

   MPStateComponent Left, Right;
};

inline
SuperblockMultiply::SuperblockMultiply(MPStateComponent const& Left_,
                                       MPStateComponent const& Right_)
   : Left(Left_), Right(Right_)
{
}

PStream::opstream& operator<<(PStream::opstream& out, DMRG const& d)
{
   return out << d.HamMatrices
              << d.Psi
              << d.Ham
              << d.HamSquared
              << d.Ident
              << d.Ortho
              << d.PsiOrthoProjector
              << d.PsiPrevC

              << d.SaveDiscardedWavefunction
              << d.LastOverlap
              << d.IsPsiConverged
              << d.IsConvergedValid
              << d.TestConverged
              << d.KeepList

              << d.TotalSweepNumber
              << d.TotalSweepRecNumber
              << d.TotalNumIterations
              << d.TotalNumMultiplies

              << d.SweepNumIterations
              << d.SweepSumStates
              << d.SweepMaxStates
              << d.SweepNumMultiplies
              << d.SweepEnergy
              << d.SweepTruncation
              << d.SweepTruncationFrac
              << d.SweepEntropy
              << d.SweepStartTime
              << d.SweepTruncatedEnergy
              << d.SweepEnergyError
              << d.SweepLastMixFactor

              << d.IterationNumMultiplies
              << d.IterationNumStates
              << d.IterationEnergy
              << d.IterationTruncation
              << d.IterationTruncationFrac
              << d.IterationEntropy
              << d.IterationEnergyBeforeTrunc
              << d.IterationEnergyVec;
}

PStream::ipstream& operator>>(PStream::ipstream& in, DMRG& d)
{
   return in >> d.HamMatrices
             >> d.Psi
             >> d.Ham
             >> d.HamSquared
             >> d.Ident
             >> d.Ortho
             >> d.PsiOrthoProjector
             >> d.PsiPrevC

             >> d.SaveDiscardedWavefunction
             >> d.LastOverlap
             >> d.IsPsiConverged
             >> d.IsConvergedValid
             >> d.TestConverged
             >> d.KeepList

             >> d.TotalSweepNumber
             >> d.TotalSweepRecNumber
             >> d.TotalNumIterations
             >> d.TotalNumMultiplies

             >> d.SweepNumIterations
             >> d.SweepSumStates
             >> d.SweepMaxStates
             >> d.SweepNumMultiplies
             >> d.SweepEnergy
             >> d.SweepTruncation
             >> d.SweepTruncationFrac
             >> d.SweepEntropy
             >> d.SweepStartTime
             >> d.SweepTruncatedEnergy
             >> d.SweepEnergyError
             >> d.SweepLastMixFactor

             >> d.IterationNumMultiplies
             >> d.IterationNumStates
             >> d.IterationEnergy
             >> d.IterationTruncation
             >> d.IterationTruncationFrac
             >> d.IterationEntropy
             >> d.IterationEnergyBeforeTrunc
             >> d.IterationEnergyVec;
}

DMRG::DMRG(CenterWavefunction const& Psi_, SplitOperator const& Ham_, bool Verbose)
   : Psi(Psi_), Ham(Ham_),
     HamSquared(prod(Ham_, Ham_, Ham_.TransformsAs())),
     Ident(Psi_.GetSymmetryList()), SaveDiscardedWavefunction(false),
     NormalizeWavefunction(false),
     IsPsiConverged(false), IsConvergedValid(false),
     MixUseEnvironment(false), UseDGKS(false), Solver("lanczos"),
     TwoStageTruncation(false), ExtraFrac(0)

{
   // construct the HamMatrix elements for the right hand side

   typedef MPStateComponent Component;

   // The initial matrix elements represent the vacuum basis on the right hand side
   BasisList Vacuum = make_vacuum_basis(Psi.GetSymmetryList());

   if (Verbose)
      std::cout << "Constructing Hamiltonian block operators..." << std::endl;

   MPStateComponent Elem(Vacuum, VectorBasis(Vacuum), VectorBasis(Vacuum));
   Elem[0](0,0) = LinearAlgebra::Matrix<double>(1,1,1);
   HamMatrices.PushRight(Elem); // add the vacuum operator, for convenience
   // apply the Right matrices
   for (int Loc = 0; Loc < Psi.RightSize(); ++Loc)
   {
      if (Verbose)
         std::cout << "site " << (Psi.RightSize()-Loc) << std::endl;
      Elem = operator_prod(Ham.LookupRight(Loc),
                              Psi.LookupRight(Loc),
                              Elem,
                              herm(Psi.LookupRight(Loc)));
      HamMatrices.PushRight(Elem);
   }

   //Now start from the left hand side and work inwards
   BasisList LeftVacuum(Psi.GetSymmetryList());
   CHECK(Psi.LookupLeft(0).Basis1().size() == 1);
   LeftVacuum.push_back(Psi.LookupLeft(0).Basis1()[0]);

   Elem = MPStateComponent(Vacuum, VectorBasis(LeftVacuum), VectorBasis(LeftVacuum));
   Elem[0](0,0) = LinearAlgebra::Matrix<double>(1,1,1);
   HamMatrices.PushLeft(Elem); // add the vacuum operator, for convenience

   for (int Loc = 0; Loc < Psi.LeftSize(); ++Loc)
   {
      Elem = operator_prod(herm(Ham.LookupLeft(Loc)),
                           herm(Psi.LookupLeft(Loc)),
                           Elem,
                           Psi.LookupLeft(Loc));
      HamMatrices.PushLeft(Elem);
   }

   // Make sure the initial wavefunction is normalized
   Psi.Normalize();

   // clear the global statistics
   TotalSweepNumber = 0;
   TotalSweepRecNumber = 0;
   TotalNumIterations = 0;
   TotalNumMultiplies = 0;

  this->DebugCheckBasis();
}

void
DMRG::AddOrthogonalState(CenterWavefunction x)
{
   // shift to line up with Psi
   while (x.LeftSize() > Psi.LeftSize())
      x.RotateLeft();

   while (x.RightSize() > Psi.RightSize())
      x.RotateRight();

   // Set up the matrix elements of the projectors
   TransformOperator Proj;

   BasisList Vacuum = make_vacuum_basis(Psi.GetSymmetryList());
   MatrixOperator VacIdent = MatrixOperator(VectorBasis(Vacuum), VectorBasis(Vacuum), Ident);
   VacIdent(0,0) = LinearAlgebra::Matrix<double>(1,1,1);
   Proj.PushRight(VacIdent);
   // apply the Right matrices
   for (int Loc = 0; Loc < x.RightSize(); ++Loc)
   {
      Proj.PushRight(operator_prod(Psi.LookupRight(Loc),
                                   Proj.Right(),
                                   herm(x.LookupRight(Loc))));
   }

   //Now start from the left hand side and work inwards
   BasisList LeftVacuum(x.GetSymmetryList());
   CHECK_EQUAL(x.LookupLeft(0).Basis1().size(), 1);
   LeftVacuum.push_back(x.LookupLeft(0).Basis1()[0]);

   VacIdent = MatrixOperator(VectorBasis(LeftVacuum), VectorBasis(LeftVacuum), Ident);
   VacIdent(0,0) = LinearAlgebra::Matrix<double>(1,1,1);

   Proj.PushLeft(VacIdent);

   for (int Loc = 0; Loc < x.LeftSize(); ++Loc)
   {
      Proj.PushLeft(operator_prod(herm(Psi.LookupLeft(Loc)),
                                  Proj.Left(),
                                  x.LookupLeft(Loc)));
   }

   Ortho.push_back(x);
   PsiOrthoProjector.push_back(Proj);

   this->DebugCheckBasis();
}

CenterWavefunction&
DMRG::Wavefunction()
{
   Psi.AttributesMutable()["Energy"] = boost::lexical_cast<std::string>(this->Energy());
   return Psi;
}

CenterWavefunction const&
DMRG::Wavefunction() const
{
   Psi.AttributesMutable()["Energy"] = boost::lexical_cast<std::string>(this->Energy());
   return Psi;
}

char ToLower(char c)
{
   return std::tolower(c);
}

void DMRG::CreateLogFiles(std::string const& BasePath, ConfList const& Conf)
{
   bool const ClearLogs = Conf.Get("ClearLogs", true);
   std::ios_base::openmode Mode = ClearLogs ? (std::ios_base::out | std::ios_base::trunc)
      : (std::ios_base::out | std::ios_base::app);

   std::string const EnergyLogFile = BasePath + ".e";
   int const EnergyLogLevel = Conf.Get("EnergyLogLevel", 0);
   Logger("EnergyLog").SetThreshold(EnergyLogLevel);
   if (EnergyLogLevel == 0 && ClearLogs)
      remove(EnergyLogFile.c_str());
   else if (EnergyLogLevel > 0)
   {
      EnergyLog = boost::shared_ptr<std::ofstream>
         (new std::ofstream(EnergyLogFile.c_str(), Mode));
      EnergyLog->precision(getenv_or_default("MP_PRECISION", 14));
      (*EnergyLog) << "#TotalSweepNumber #SweepNumber #LeftSize #RightSize "
                   << "#NumStates #NumMulti #Trunc #Trunc50 #Entropy #Energy\n";
      Logger("EnergyLog").SetStream(*EnergyLog);
   }

   std::string const SweepLogFile = BasePath + ".sweep";
   int const SweepLogLevel = Conf.Get("SweepLogLevel", 0);
   Logger("SweepLog").SetThreshold(SweepLogLevel);
   if (SweepLogLevel == 0 && ClearLogs)
      remove(SweepLogFile.c_str());
   else if (SweepLogLevel > 0)
   {
      SweepLog = boost::shared_ptr<std::ofstream>
         (new std::ofstream(SweepLogFile.c_str(), Mode));
      SweepLog->precision(getenv_or_default("MP_PRECISION", 14));
      (*SweepLog) << "#TotalSweepNum #SweepNum #NIterations #States "
         "#NMult #WTime #Trunc #E_trunc #E_sigma #Entropy #(H-E)^2 #Overlap #ODiff #Energy #MixFactor #Converged\n";
      Logger("SweepLog").SetStream(*SweepLog);
   }

   std::string const DensityLogFile = BasePath + ".density";
   int const DensityLogLevel = Conf.Get("DensityLogLevel", 0);
   Logger("DensityLog").SetThreshold(DensityLogLevel);
   if (DensityLogLevel == 0 && ClearLogs)
      remove(DensityLogFile.c_str());
   else if (DensityLogLevel > 0)
   {
      DensityLog = boost::shared_ptr<std::ofstream>
         (new std::ofstream(DensityLogFile.c_str(), Mode));
      DensityLog->precision(getenv_or_default("MP_PRECISION", 14));
      Logger("DensityLog").SetStream(*DensityLog);
   }

   // CpuLog
   // DiagLog

   ConvergenceOverlapTruncationScale = Conf.Get("Convergence::OverlapTruncationRatio", 1.0);
   ConvergenceOverlapDifferenceOverlapScale = Conf.Get("Convergence::OverlapDerivativeRatio", 1.0);
   ConvergenceSweepTruncMin = Conf.Get("Convergence::TruncationCutoff", 1E-13);
   ConvergenceOverlapMin = Conf.Get("Convergence::OverlapCutoff", 1E-13);
   MixUseEnvironment = Conf.Get("MixUseEnvironment", false);
   SaveDiscardedWavefunction = Conf.Get("SaveDiscardedWavefunction", false);
   NormalizeWavefunction = Conf.Get("Normalize", true);
   Solver = Conf.Get("Solver", "Lanczos");
   std::transform(Solver.begin(), Solver.end(), Solver.begin(), &ToLower);
   TwoStageTruncation = Conf.Get("TwoStageTruncation", false);
   ExtraFrac = Conf.Get("ExtraFrac", 0.2);
   UseDGKS = Conf.Get("UseDGKS", false);
}

void DMRG::RestoreLogFiles(std::string const& BasePath, ConfList const& Conf)
{
   std::string const EnergyLogFile = BasePath + ".e";
   int const EnergyLogLevel = Conf.Get("EnergyLogLevel", 0);
   Logger("EnergyLog").SetThreshold(EnergyLogLevel);
   if (EnergyLogLevel > 0)
   {
      EnergyLog = boost::shared_ptr<std::ofstream>
         (new std::ofstream(EnergyLogFile.c_str(), std::ios_base::out | std::ios_base::app));
      EnergyLog->precision(getenv_or_default("MP_PRECISION", 14));
      Logger("EnergyLog").SetStream(*EnergyLog);
   }

   std::string const SweepLogFile = BasePath + ".sweep";
   int const SweepLogLevel = Conf.Get("SweepLogLevel", 0);
   Logger("SweepLog").SetThreshold(SweepLogLevel);
   if (SweepLogLevel > 0)
   {
      SweepLog = boost::shared_ptr<std::ofstream>
         (new std::ofstream(SweepLogFile.c_str(), std::ios_base::out | std::ios_base::app));
      SweepLog->precision(getenv_or_default("MP_PRECISION", 14));
      Logger("SweepLog").SetStream(*SweepLog);
   }

   std::string const DensityLogFile = BasePath + ".density";
   int const DensityLogLevel = Conf.Get("DensityLogLevel", 0);
   Logger("DensityLog").SetThreshold(DensityLogLevel);
   if (DensityLogLevel > 0)
   {
      DensityLog = boost::shared_ptr<std::ofstream>
         (new std::ofstream(DensityLogFile.c_str(), std::ios_base::out | std::ios_base::app));
      DensityLog->precision(getenv_or_default("MP_PRECISION", 14));
      Logger("DensityLog").SetStream(*DensityLog);
   }

   // CpuLog
   // DiagLog

   ConvergenceOverlapTruncationScale = Conf.Get("Convergence::OverlapTruncationRatio", 1.0);
   ConvergenceOverlapDifferenceOverlapScale = Conf.Get("Convergence::OverlapDerivativeRatio", 1.0);
   ConvergenceSweepTruncMin = Conf.Get("Convergence::TruncationCutoff", 1E-13);
   ConvergenceOverlapMin = Conf.Get("Convergence::OverlapCutoff", 1E-13);
   MixUseEnvironment = Conf.Get("MixUseEnvironment", false);
   SaveDiscardedWavefunction = Conf.Get("SaveDiscardedWavefunction", false);
   NormalizeWavefunction = Conf.Get("Normalize", true);
   Solver = Conf.Get("Solver", "Lanczos");
   std::transform(Solver.begin(), Solver.end(), Solver.begin(), &ToLower);
   TwoStageTruncation = Conf.Get("TwoStageTruncation", false);
   ExtraFrac = Conf.Get("ExtraFrac", 0.0);
   UseDGKS = Conf.Get("UseDGKS", false);

   // debugging; set the precision for cerr
   std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

   std::cout << "Convergence::OverlapTruncationRatio = " << ConvergenceOverlapTruncationScale << '\n';
   std::cout << "Convergence::OverlapDerivativeRatio = " << ConvergenceOverlapDifferenceOverlapScale << '\n';
   std::cout << "Convergence::TruncationCutoff = " << ConvergenceSweepTruncMin << '\n';
   std::cout << "Convergence::OverlapCutoff = " << ConvergenceOverlapMin << '\n';
   std::cout << "MixUseEnvironment = " << MixUseEnvironment << '\n';
   std::cout << "SaveDiscardedWavefunction = " << SaveDiscardedWavefunction << '\n';
   std::cout << "Normalize = " << NormalizeWavefunction << '\n';
   std::cout << "Solver = " << Solver << '\n';
   std::cout << "TwoStageTruncation = " << TwoStageTruncation << '\n';
   std::cout << "ExtraFrac = " << ExtraFrac << '\n';
   std::cout << "UseDGKS = " << UseDGKS << '\n';
}

void DMRG::StartIteration()
{
   IterationNumMultiplies = 0;
   IterationNumStates = 0;
   IterationEnergy = 0;
   IterationTruncation = 0;
   IterationTruncationFrac = 0;
   IterationEntropy = 0;
}

void DMRG::EndIteration()
{
   // recalculate the energy, as it will have changed after truncation
   if (IterationTruncation != 0)
      IterationEnergy = this->Energy();
   msg_log(1, "EnergyLog") << TotalSweepNumber << ' '
                           << TotalSweepRecNumber << ' '
                           << this->LeftSize() << ' '
                           << this->RightSize() << ' '
                           << IterationNumStates << ' '
                           << IterationNumMultiplies << ' '
                           << IterationTruncation << ' '
                           << IterationTruncationFrac << ' '
                           << IterationEntropy << ' '
                           << IterationEnergy << '\n';

   ++SweepNumIterations;
   SweepNumMultiplies += IterationNumMultiplies;
   SweepSumStates += IterationNumStates;
   SweepMaxStates = std::max(SweepMaxStates, IterationNumStates);
   SweepEnergy = std::min(SweepEnergy, IterationEnergy);
   SweepTruncatedEnergy += IterationEnergyBeforeTrunc - IterationEnergy;
   SweepTruncation += IterationTruncation;
   SweepTruncationFrac += IterationTruncationFrac;
   SweepEntropy = std::max(SweepEntropy, IterationEntropy);
   IterationEnergyVec.push_back(IterationEnergy);
}

void DMRG::StartSweep(bool IncrementSweepNumber, double /* Broad_ */)
{
   ++TotalSweepNumber;
   if (IncrementSweepNumber)
      ++TotalSweepRecNumber;

   SweepNumIterations = 0;
   SweepSumStates = 0;
   SweepMaxStates = 0;
   SweepNumMultiplies = 0;
   SweepEnergy = 1E100;
   SweepTruncation = 0;
   SweepTruncationFrac = 0;
   SweepEntropy = 0;
   SweepStartTime = ProcControl::GetCumulativeElapsedTime();
   SweepTruncatedEnergy = 0;
   TestConverged = false;
   IterationEnergyVec.clear();

   Psi.Normalize();

   PsiPrevC = Psi.Center();

   // Initialize the keep list
   KeepList.clear();

   IsConvergedValid = false;
}

double DMRG::FidelityLoss() const
{
   return 2.0 * (1.0 - norm_frob(inner_prod(Psi.Center(), PsiPrevC))
                 / norm_frob(Psi.Center()));
}

void DMRG::EndSweep()
{
   TotalNumIterations += SweepNumIterations;
   TotalNumMultiplies += SweepNumMultiplies;

   SweepEnergyError = std::sqrt(statistics::variance(IterationEnergyVec.begin(),
                                                     IterationEnergyVec.end())
                                / SweepNumIterations);

   double HE2 = -1;

   double Overlap = 2.0 * (1.0 - norm_frob(inner_prod(Psi.Center(), PsiPrevC))
                           / norm_frob(Psi.Center()));
   double OverlapDifference = 2;
   if (LastOverlap)
   {
      OverlapDifference = std::abs(*LastOverlap - Overlap);
      IsPsiConverged = ((Overlap < ConvergenceOverlapTruncationScale * SweepTruncation)
                        || SweepTruncation < ConvergenceSweepTruncMin)
         && ((OverlapDifference < ConvergenceOverlapDifferenceOverlapScale * Overlap)
             || Overlap < ConvergenceOverlapMin);
      IsConvergedValid = true;
   }

   LastOverlap = Overlap;


   // We only calculate HE2 if we are converged
   if (this->IsConverged() || TestConverged)
   {
      double E2 = expectation(Psi, HamSquared, Psi).real();
      double E = this->Energy();
      HE2 = (E2 - E*E) / norm_frob_sq(Psi.Center());
   }
   else
   {
      HE2 = 0;
   }

   int Converged = this->IsConverged();
   double SweepTime = ProcControl::GetCumulativeElapsedTime() - SweepStartTime;

   msg_log(1, "SweepLog") << TotalSweepNumber << ' '
                          << TotalSweepRecNumber << ' '
                          << SweepNumIterations << ' '
                          << SweepMaxStates << ' '
                          << SweepNumMultiplies << ' '
                          << SweepTime << ' '
                          << std::max(SweepTruncation, 0.0) << ' '
                          << SweepTruncatedEnergy << ' '
                          << SweepEnergyError << ' '
                          << SweepEntropy << ' '
                          << HE2 << ' '
                          << Overlap << ' '
                          << OverlapDifference << ' '
                          << SweepEnergy << ' '
                          << SweepLastMixFactor << ' '
                          << Converged
                          << std::endl;
   msg_log(1, "EnergyLog") << std::flush;
}

double DMRG::Energy() const
{
   //   DEBUG_TRACE(HamMatrices.Left())(HamMatrices.Right());
   //   DEBUG_TRACE(scalar_prod(Psi.Center(), Psi.Center()));
   //   DEBUG_TRACE(scalar_prod(PsiP, PsiP));
   return inner_prod(operator_prod(HamMatrices.Left(),
                                   Psi.Center(),
                                   herm(HamMatrices.Right())),
                     Psi.Center()).real();
}

int OperatorMaxDimension(MatrixOperator const& x)
{
   int Dim = 0;
   for (unsigned i = 0; i < x.Basis1().size(); ++i)
   {
      for (unsigned j = 0; j < x.Basis2().size(); ++j)
      {
         if (is_transform_target(x.Basis2()[j], x.TransformsAs(), x.Basis1()[i]))
            Dim += x.Basis1().dim(i) * x.Basis2().dim(j);
      }
   }
   return Dim;
}

struct Preconditioner
{
   typedef std::complex<double> result_type;
   typedef std::complex<double> first_argument_type;
   typedef std::complex<double> second_argument_type;

   Preconditioner(double Theta_) : Theta(Theta_) {}

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      return x / (y - Theta);
   }

   double Theta;
};

void Precondition(MatrixOperator& r, MatrixOperator const& Diag, double Theta)
{
   for (LinearAlgebra::iterator<MatrixOperator>::type I = iterate(r); I; ++I)
   {
      for (LinearAlgebra::inner_iterator<MatrixOperator>::type J = iterate(I); J; ++J)
      {
         LinearAlgebra::const_inner_iterator<MatrixOperator>::type DiagJ =
            iterate_at(Diag.data(), J.index1(), J.index2());
         *J = LinearAlgebra::transform(*J, *DiagJ, Preconditioner(Theta));
      }
   }
}

double DMRG::Solve(int MaxIterations)
{
   // 31-12-2006: Fix the sign of the resulting eigenvector.
   // TODO: What happens in the complex case?
   DEBUG_TRACE("DMRG::Solve");
   MatrixOperator PsiOld = Psi.Center();

   // if the wavefunction is zero, add a random component
   double pNorm = norm_frob_sq(Psi.Center());
   if (pNorm < 1E-12)
   {
      Psi.Center() += MakeRandomMatrixOperator(Psi.Center().Basis1(), Psi.Center().Basis2());
      Psi.Center() *= 1.0 / norm_frob(Psi.Center());
   }

   //   TRACE(Psi.Center())(norm_frob(Psi.Center()));
   //   DEBUG_TRACE(Psi.Center().Basis1())(Psi.Center().Basis2());
   int Iterations = std::max(1, std::min(MaxIterations,
                                         OperatorMaxDimension(Psi.Center()) - int(Ortho.size())));

   std::vector<MatrixOperator> OrthoSet(Ortho.size());
   for (std::size_t j = 0; j < OrthoSet.size(); ++j)
   {
      OrthoSet[j] = triple_prod(PsiOrthoProjector[j].Left(),
                                Ortho[j].Center(),
                                herm(PsiOrthoProjector[j].Right()));
   }

   double Energy=0;
   if (Solver == "lanczos")
   {
      Energy = Lanczos(Psi.Center(),
                       SuperblockMultiply(HamMatrices.Left(),
                                          HamMatrices.Right()),
                       Iterations,
                       OrthoSet,
                       UseDGKS);
   }
   else if (Solver == "davidson")
   {
      CHECK(OrthoSet.empty())("Davidson solver does not support orthogonal states yet.");
      MatrixOperator Diag = extract_diagonal(HamMatrices.Left(), herm(HamMatrices.Right()));
      Energy =  LinearSolvers::Davidson(Psi.Center(),
                                        Diag,
                                        SuperblockMultiply(HamMatrices.Left(),
                                                           HamMatrices.Right()),
                                        Iterations);
   }
   else if (Solver == "arnoldi")
   {
      CHECK(OrthoSet.empty())("Arnoldi solver does not support orthogonal states yet.");
      std::complex<double> Eigen = LinearSolvers::Arnoldi(Psi.Center(),
                                                          SuperblockMultiply(HamMatrices.Left(),
                                                                             HamMatrices.Right()),
                                                          Iterations);
      Energy = Eigen.real();
      //      EnergyImag = Eigen.imag();
   }
   else
      PANIC("Unknown solver")(Solver);

   IterationEnergy = Energy;
   IterationNumMultiplies = Iterations;

   // sign of the eigenvector
   if (inner_prod(PsiOld, Psi.Center()).real() < 0)
      Psi.Center() *= -1.0;

   //CHECK_CLOSE(Energy, this->Energy());

   IterationEnergyBeforeTrunc = Energy;

   return Energy;
}

void DMRG::PrepareConvergenceTest()
{
   TestConverged = true;
}

bool DMRG::IsConverged() const
{
   return IsConvergedValid && IsPsiConverged;
}

void DMRG::ExpandLeft()
{
   DEBUG_TRACE("DMRG::ExpandLeft");
   MPStateComponent PsiLeftOld = Psi.Left();
   Psi.Left() = prod(Psi.Left(), Psi.Center());
   Psi.Center() = ExpandBasis2(Psi.Left());

   MatrixOperator U = scalar_prod(herm(Psi.Left()), PsiLeftOld);
   PsiPrevC = U * PsiPrevC;

   HamMatrices.PopLeft();
   MPStateComponent New = operator_prod(herm(Ham.Left()),
                                        herm(Psi.Left()),
                                        HamMatrices.Left(),
                                        (Psi.Left()));
   HamMatrices.PushLeft(New);

   // Expand the projectors for the othogonal set
   for (std::size_t j = 0; j < Ortho.size(); ++j)
   {
      PsiOrthoProjector[j].PopLeft();
      PsiOrthoProjector[j].PushLeft(operator_prod(herm(Psi.Left()),
                                               PsiOrthoProjector[j].Left(),
                                               Ortho[j].Left()));
   }
  this->DebugCheckBasis();
}

void DMRG::ExpandRight()
{
   DEBUG_TRACE("DMRG::ExpandRight");
   MPStateComponent PsiRightOld = Psi.Right();
   Psi.Right() = prod(Psi.Center(), Psi.Right());
   Psi.Center() = ExpandBasis1(Psi.Right());

   MatrixOperator U = scalar_prod(PsiRightOld, herm(Psi.Right()));
   PsiPrevC = PsiPrevC * U;

   HamMatrices.PopRight();
   MPStateComponent New = operator_prod(Ham.Right(),
                                        Psi.Right(),
                                        HamMatrices.Right(),
                                        herm(Psi.Right()));
   HamMatrices.PushRight(New);

   // Expand the projectors for the othogonal set
   for (std::size_t j = 0; j < Ortho.size(); ++j)
   {
      PsiOrthoProjector[j].PopRight();
      PsiOrthoProjector[j].PushRight(operator_prod(Psi.Right(),
                                                PsiOrthoProjector[j].Right(),
                                                herm(Ortho[j].Right())));
   }
  this->DebugCheckBasis();
}

void DMRG::ShiftRightAndExpand()
{
   DEBUG_TRACE("DMRG::ShiftRightAndExpand");
   MPStateComponent OldPsiRight = Psi.Right();
   Psi.PushLeft(prod(Psi.Center(), Psi.Right()));
   Psi.PopRight();

   // Expand the new left basis, this is the 'old DMRG' equivalent
   // of inverting the truncation operator at the end of
   // the wavefunction transformation
   Psi.Center() = ExpandBasis2(Psi.Left());

   // rotate the Ham operator

   //   TRACE(Ham.Center())(Ham.Left())(Ham.Right());

   PsiPrevC = operator_prod(herm(Psi.Left()), PsiPrevC, OldPsiRight);

   Ham.RotateRight();

   //   TRACE(Ham.Center())(Ham.Left())(Ham.Right());

   // calculate the new Ham matrix and rotate
   MPStateComponent New = operator_prod(herm(Ham.Left()),
                                        herm(Psi.Left()),
                                        HamMatrices.Left(),
                                        Psi.Left());
   HamMatrices.PushLeft(New);
   HamMatrices.PopRight();

   // Expand the projectors for the othogonal set
   for (std::size_t j = 0; j < Ortho.size(); ++j)
   {
      Ortho[j].RotateRight();
      PsiOrthoProjector[j].PushLeft(operator_prod(herm(Psi.Left()),
                                                  PsiOrthoProjector[j].Left(),
                                                  Ortho[j].Left()));
      PsiOrthoProjector[j].PopRight();
   }
  this->DebugCheckBasis();
}

void DMRG::ShiftLeftAndExpand()
{
   DEBUG_TRACE("DMRG::ShiftLeftAndExpand");
   MPStateComponent OldPsiLeft = Psi.Left();
   Psi.PushRight(prod(Psi.Left(), Psi.Center()));
   Psi.PopLeft();

   // 'invert' the truncation on the new right block
   Psi.Center() = ExpandBasis1(Psi.Right());

   PsiPrevC = operator_prod(OldPsiLeft, PsiPrevC, herm(Psi.Right()));

   // rotate the Ham operator

   //   TRACE(Ham.Center())(Ham.Left())(Ham.Right());


   Ham.RotateLeft();

   //   TRACE(Ham.Center())(Ham.Left())(Ham.Right());

   // calculate the new matrices
   MPStateComponent New = operator_prod(Ham.Right(),
                                           Psi.Right(),
                                           HamMatrices.Right(),
                                           herm(Psi.Right()));
   HamMatrices.PushRight(New);
   HamMatrices.PopLeft();

   // Expand the projectors for the othogonal set
   for (std::size_t j = 0; j < Ortho.size(); ++j)
   {
      Ortho[j].RotateLeft();
      PsiOrthoProjector[j].PushRight(operator_prod(Psi.Right(),
                                                PsiOrthoProjector[j].Right(),
                                                herm(Ortho[j].Right())));
      PsiOrthoProjector[j].PopLeft();
   }
  this->DebugCheckBasis();
}

// Two-stage truncation scheme.  Rho is the ordinary density matrix, that
// is used in the first round truncation according to TruncationSpec.
// Rho2 is used within the subspace of discarded states for a second round,
// where the kept states are chosen according to the desired 'expansion fraction',
// from Spec.EXtraFrac().
MatrixOperator TruncateTwoStage(MatrixOperator const& Rho, MatrixOperator const& Rho2,
                                StatesInfo const& SInfo,
                                double MixFactor,
                                double ExtraFrac,
                                KeepListType& KeepList,
                                std::set<QuantumNumbers::QuantumNumber> const& AddedQN,
                                TruncationInfo& Info)
{
   // First round: construct the density matrix
   DensityMatrix<MatrixOperator> DM(Rho);
   DensityMatrix<MatrixOperator>::const_iterator DMPivot =
      TruncateFixTruncationErrorAbsolute(DM.begin(),
                                         DM.end(),
                                         SInfo,
                                         Info);
   // Split into kept and discarded states
   std::list<EigenInfo> KeptStates(DM.begin(), DMPivot);
   std::list<EigenInfo> DiscardStates(DMPivot, DM.end());
   // Update the keep list.  It would perhaps be better to do this with respect
   // to the stage 2 density matrix, but easier to do it here
   UpdateKeepList(KeepList,
                  AddedQN,
                  DM.Basis(),
                  KeptStates,
                  DiscardStates,
                  Info);

   // Dimensions of each quantum number subspace
   EigenDimensionsType EDim = EigenDimensions(Rho.Basis1(), KeptStates.begin(), KeptStates.end());
   // Get the desired size of additional fraction per subspace
   UpdateEigenDimensions(EDim, ExtraFrac);

   // Projector onto the kept and discarded states
   MatrixOperator UDiscard = DM.ConstructTruncator(DiscardStates.begin(), DiscardStates.end());
   MatrixOperator UKeep = DM.ConstructTruncator(KeptStates.begin(), KeptStates.end());
   // Second round
   // Project Rho2 onto the discarded states
   MatrixOperator Rho2Trunc = triple_prod(UDiscard, Rho2, herm(UDiscard));

   DensityMatrix<MatrixOperator> DM2(Rho2Trunc);

   //TRACE("DM");
   //DM.DensityMatrixReport(std::cout);
   //TRACE("DM2");
   //DM2.DensityMatrixReport(std::cout);

   TruncationInfo Info2;
   std::list<EigenInfo> SecondStage = TruncateFixEigenDimensions(EDim, Rho2Trunc.Basis1(),
                                                                 int(SInfo.MinStates * ExtraFrac+0.5),
                                                                 DM2.begin(), DM2.end(),
                                                                 Info2);

   Info.KeptStates_ += Info2.KeptStates_;

   // Second stage projector
   MatrixOperator USecondStage = DM2.ConstructTruncator(SecondStage.begin(), SecondStage.end());
   // convert back to an operator over the complete basis
   USecondStage = USecondStage * UDiscard;

   // Now we construct the tensor sum of the first stage and second stage kept states
   MatrixOperator UCombined = tensor_col_sum(UKeep, USecondStage);

   return UCombined;
}

TruncationInfo DMRG::TruncateLeft(StatesInfo const& SInfo, double CFactor)
{
   DEBUG_TRACE("DMRG::TruncateLeft")(SInfo);
   MatrixOperator Rho = scalar_prod(Psi.Center(), herm(Psi.Center()));

   MatrixOperator U; // the truncation operator
   TruncationInfo Info;

   SweepLastMixFactor = CFactor;  // so we can put the mix factor into the .sweep file

   if (TwoStageTruncation)
   {
      MatrixOperator Rho2 = Rho;
      if (CFactor != 0)
      {
         for (unsigned i = 0; i < HamMatrices.Left().size(); ++i)
         {
            MatrixOperator Correction = triple_prod(HamMatrices.Left()[i],
                                                    Rho,
                                                    herm(HamMatrices.Left()[i]));
            std::complex<double> Tr = trace(Correction);
            if (norm_2(Tr) > std::numeric_limits<double>::epsilon() * 100)
            {
               Correction *= CFactor / (Tr * double(HamMatrices.Left().size()));
               Rho2 += Correction;
            }
         }
      }

      // Truncation specification, get ExtraFrac from *this
      U = TruncateTwoStage(Rho, Rho2, SInfo, CFactor, ExtraFrac, KeepList,
                           QuantumNumbersInBasis(adjoint(Psi.Left().SiteBasis())),
                           Info);
   }
   else // old truncation scheme
   {
      if (CFactor != 0)
      {
         MatrixOperator Correction;
         if (MixUseEnvironment)
         {
            MatrixOperator RhoR = scalar_prod(herm(Psi.Center()), Psi.Center());
            Correction =
               operator_prod(trace_prod(herm(HamMatrices.Right()),
                                        prod(RhoR, HamMatrices.Right())),
                             HamMatrices.Left(),
                             Rho,
                             herm(HamMatrices.Left()));
         }
         else
         {
            for (unsigned i = 0; i < HamMatrices.Left().size(); ++i)
            {
               MatrixOperator Correction = triple_prod(HamMatrices.Left()[i],
                                                       Rho,
                                                       herm(HamMatrices.Left()[i]));
               std::complex<double> Tr = trace(Correction);
               if (norm_2(Tr) > std::numeric_limits<double>::epsilon() * 100)
               {
                  Correction *= (CFactor / trace(Correction) * double(HamMatrices.Left().size()));
                  Rho += Correction;
               }
            }
         }
      }

      DensityMatrix<MatrixOperator> DM(Rho);
      DensityMatrix<MatrixOperator>::const_iterator DMPivot =
         TruncateFixTruncationErrorAbsolute(DM.begin(), DM.end(),
                                            SInfo,
                                            Info);

      std::list<EigenInfo> KeptStates(DM.begin(), DMPivot);
      std::list<EigenInfo> DiscardStates(DMPivot, DM.end());
      UpdateKeepList(KeepList,
                     QuantumNumbersInBasis(adjoint(Psi.Left().SiteBasis())),
                     DM.Basis(),
                     KeptStates,
                     DiscardStates,
                     Info);

      U = DM.ConstructTruncator(KeptStates.begin(), KeptStates.end());
   }

   IterationNumStates = Info.KeptStates();
   IterationTruncation += Info.TruncationError();
   IterationEntropy = std::max(IterationEntropy, Info.KeptEntropy());

   // truncate the wavefunction
   Psi.Center() = U * Psi.Center();
   Psi.Left() = prod(Psi.Left(), herm(U));
   Psi.Normalize();

   PsiPrevC = U*PsiPrevC;

   // truncate the Ham matrix
   HamMatrices.Left() = triple_prod(U, HamMatrices.Left(), herm(U));

   // truncate the ortho projections
   for (std::size_t j = 0; j < Ortho.size(); ++j)
   {
      PsiOrthoProjector[j].Left() = prod(U, PsiOrthoProjector[j].Left(), Ident);
   }

   this->DebugCheckBasis();
   return Info;
}

TruncationInfo DMRG::TruncateRight(StatesInfo const& SInfo, double CFactor)
{
   DEBUG_TRACE("DMRG::TruncateRight")(SInfo)(CFactor);
   MatrixOperator Rho = scalar_prod(herm(Psi.Center()), Psi.Center());
   MatrixOperator U; // the truncation operator
   TruncationInfo Info;

   SweepLastMixFactor = CFactor;  // so we can put the mix factor into the .sweep file

   if (TwoStageTruncation)
   {
      MatrixOperator Rho2 = Rho;
      if (CFactor != 0)
      {
         for (unsigned i = 0; i < HamMatrices.Right().size(); ++i)
         {
            MatrixOperator Correction = triple_prod(HamMatrices.Right()[i],
                                                    Rho,
                                                    herm(HamMatrices.Right()[i]));
            std::complex<double> Tr = trace(Correction);
            if (norm_2(Tr) > std::numeric_limits<double>::epsilon() * 100)
            {
               Correction *= CFactor / trace(Correction) * double(HamMatrices.Right().size());
               Rho2 += Correction;
            }
         }
      }

      // Truncation specification, get ExtraFrac from *this
      U = TruncateTwoStage(Rho, Rho2, SInfo, CFactor, ExtraFrac, KeepList,
                           QuantumNumbersInBasis(Psi.Right().SiteBasis()),
                           Info);
   }
   else // old truncation scheme
   {
      if (CFactor != 0)
      {
         MatrixOperator Correction;
         if (MixUseEnvironment)
         {
            MatrixOperator RhoL = scalar_prod(Psi.Center(), herm(Psi.Center()));
            Correction =
               operator_prod(trace_prod(herm(HamMatrices.Left()),
                                        prod(RhoL, HamMatrices.Left())),
                             HamMatrices.Right(),
                             Rho,
                             herm(HamMatrices.Right()));
            Correction *= (CFactor / trace(Correction));
            //Correction *= 0.25;   // this fixes up the scaling problem!
            Rho += Correction;
         }
         else
         {
            for (unsigned i = 0; i < HamMatrices.Right().size(); ++i)
            {
               MatrixOperator Correction = triple_prod(HamMatrices.Right()[i],
                                                       Rho,
                                                       herm(HamMatrices.Right()[i]));
               std::complex<double> Tr = trace(Correction);
               if (norm_2(Tr) > std::numeric_limits<double>::epsilon() * 100)
               {
                  Correction *= CFactor / trace(Correction) * double(HamMatrices.Right().size());
                  Rho += Correction;
               }
            }
         }
      }

      DensityMatrix<MatrixOperator> DM(Rho);
      DensityMatrix<MatrixOperator>::const_iterator DMPivot =
         TruncateFixTruncationErrorAbsolute(DM.begin(), DM.end(),
                                            SInfo,
                                            Info);

      std::list<EigenInfo> KeptStates(DM.begin(), DMPivot);
      std::list<EigenInfo> DiscardStates(DMPivot, DM.end());
      UpdateKeepList(KeepList,
                     QuantumNumbersInBasis(Psi.Right().SiteBasis()),
                     DM.Basis(),
                     KeptStates,
                     DiscardStates,
                     Info);

      U = DM.ConstructTruncator(KeptStates.begin(), KeptStates.end());
   }

   IterationNumStates = Info.KeptStates();
   IterationTruncation += Info.TruncationError();
   IterationEntropy = std::max(IterationEntropy, Info.KeptEntropy());

#if 0
   // Discarded wavefunction
   if (SaveDiscardedWavefunction)
   {
      CenterWavefunction Discard = Psi;
      MatrixOperator UDiscard = DM.ConstructTruncator(Piv, DM.end());
      Discard.Center() = prod(Discard.Center(), adjoint(UDiscard), Ident);
      Discard.Right() = prod(UDiscard, Discard.Right());
      Discard.RotateToNormalForm();
      std::string Filename = "discard-wavefunction."
         + boost::lexical_cast<std::string>(TotalSweepNumber) + "."
         + boost::lexical_cast<std::string>(Psi.LeftSize());
      pvalue_ptr<CenterWavefunction> OutPsi(new CenterWavefunction(Discard));
      pheap::ExportHeap(Filename, OutPsi);
   }
#endif

   // truncate the wavefunction
   Psi.Center() = Psi.Center() * herm(U);
   Psi.Right() = prod(U, Psi.Right());
   if (NormalizeWavefunction)
      Psi.Normalize();

   PsiPrevC = PsiPrevC*herm(U);

   //std::cout << 1.0-trace(scalar_prod(herm(Psi.Center()), Psi.Center())).real() << '\n';

   // truncate the Ham matrix
   HamMatrices.Right() = triple_prod(U, HamMatrices.Right(), herm(U));

   // truncate the ortho projections
   for (std::size_t j = 0; j < Ortho.size(); ++j)
   {
      PsiOrthoProjector[j].Right() = prod(U, PsiOrthoProjector[j].Right(), Ident);
   }

   this->DebugCheckBasis();
   return Info;
}

#if 0
void DMRG::InsertSitesDeltaShift(std::vector<SiteBasis> const& LeftSites,
                                 std::vector<SiteBasis> const& RightSites,
                                 QuantumNumber NewTarget,
                                 SplitOperator const& NewHam)
{
   // We can't handle orthogonal states here
   CHECK(Ortho.empty());

   QuantumNumber CurrentTarget = Psi.LeftVacuumBasis()[0];
   QuantumNumbers::Projection Shift = difference(NewTarget, CurrentTarget);
   QuantumNumber q = heighest_weight(Shift);

   // Firstly, apply the delta shift to the left matrices
   DEBUG_TRACE("DeltaShift by")(q)(Shift);
   for (int i = 0; i < Psi.LeftSize(); ++i)
   {
      TRACE(i)(Psi.LookupLeft(i).Basis1())(Psi.LookupLeft(i).Basis2());
      TRACE(i)(Psi.LookupLeft(i)[0].Basis1())(Psi.LookupLeft(i)[0].Basis2());
      Psi.LookupLeft(i).delta_shift(q, Shift);
   }

   // Add the left sites
   for (unsigned i = 0; i < LeftSites.size(); ++i)
   {
      Psi.PushLeft(MPStateComponent::ConstructFullBasis2(Psi.Left().Basis2(),
                                                         LeftSites[i].Basis()));
   }

   // Add the right sites
   for (int i = int(RightSites.size())-1; i >= 0; --i)
   {
      Psi.PushRight(MPStateComponent::ConstructFullBasis1(RightSites[i].Basis(),
                                                          Psi.Right().Basis1()));
   }

   // Replace the Center matrix with some random state, since we don't have a better
   // guess.
   Psi.Center() = MakeRandomMatrixOperator(Psi.Left().Basis2(), Psi.Right().Basis1());
   Psi.Center() *= 1.0 / norm_frob(Psi.Center());  // normalize ti

   // Now we need to reconstruct the Hamiltonian matrix elements
   Ham = NewHam;
   HamSquared = prod(Ham, Ham, Ham.TransformsAs());
   InitializeSuperblockStack(HamMatrices, Psi, Ham, Psi);
}
#endif

void DMRG::InsertSites(std::vector<SiteBasis> const& LeftSites,
                       std::vector<SiteBasis> const& RightSites,
                       QuantumNumber NewTarget,
                       SplitOperator const& NewHam)
{
   // We can't handle orthogonal states here
   CHECK(Ortho.empty());

#if 0
   QuantumNumber CurrentTarget = Psi.LeftVacuumBasis()[0];
   QuantumNumbers::Projection Shift = difference(NewTarget, CurrentTarget);
   QuantumNumber q = heighest_weight(Shift);

   // Firstly, apply the delta shift to the left matrices
   for (int i = 0; i < Psi.LeftSize(); ++i)
   {
      Psi.LookupLeft(i).delta_shift(q, Shift);
   }
#endif

   // Add the left sites
   for (unsigned i = 0; i < LeftSites.size(); ++i)
   {
      Psi.PushLeft(MPStateComponent::ConstructFullBasis2(Psi.Left().Basis2(),
                                                         LeftSites[i].Basis()));
   }

   // Add the right sites
   for (int i = int(RightSites.size())-1; i >= 0; --i)
   {
      Psi.PushRight(MPStateComponent::ConstructFullBasis1(RightSites[i].Basis(),
                                                          Psi.Right().Basis1()));
   }

   // Replace the Center matrix with some random state, since we don't have a better
   // guess.
   QuantumNumber CurrentTarget = Psi.LeftVacuumBasis()[0];
   QuantumNumbers::Projection Shift = difference(NewTarget, CurrentTarget);
   QuantumNumber q = heighest_weight(Shift);

   Psi.Center() = MakeRandomMatrixOperator(Psi.Left().Basis2(), Psi.Right().Basis1(), q);
   Psi.Center() *= 1.0 / norm_frob(Psi.Center());  // normalize ti

   // Now we need to reconstruct the Hamiltonian matrix elements
   Ham = NewHam;
   HamSquared = prod(Ham, Ham, Ham.TransformsAs());
   InitializeSuperblockStack(HamMatrices, Psi, Ham, Psi);
}

void DMRG::DebugCheckBasis() const
{
   DEBUG_CHECK_EQUAL(HamMatrices.Left().Basis1(), Psi.Center().Basis1());
   DEBUG_CHECK_EQUAL(HamMatrices.Right().Basis1(), Psi.Center().Basis2());

   DEBUG_CHECK_EQUAL(HamMatrices.Left().Basis2(), Psi.Center().Basis1());
   DEBUG_CHECK_EQUAL(HamMatrices.Right().Basis2(), Psi.Center().Basis2());

   for (std::size_t j = 0; j < Ortho.size(); ++j)
   {
      DEBUG_CHECK_EQUAL(PsiOrthoProjector[j].Left().Basis1(), Psi.Center().Basis1());
      DEBUG_CHECK_EQUAL(PsiOrthoProjector[j].Right().Basis1(), Psi.Center().Basis2());

      DEBUG_CHECK_EQUAL(PsiOrthoProjector[j].Left().Basis2(), Ortho[j].Center().Basis1());
      DEBUG_CHECK_EQUAL(PsiOrthoProjector[j].Right().Basis2(), Ortho[j].Center().Basis2());
   }
}
