// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/dmrg.cpp
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
// $Id$

#include "dmrg.h"
#include "lanczos-ortho.h"
#include "davidson.h"
#include "arnoldi.h"
#include "mps/density.h"
#include "mps/truncation.h"
#include "pstream/optional.h"
#include <boost/optional.hpp>
#include <boost/none.hpp>
#include "common/statistics.h"
#include "common/environment.h"
#include "tensor/tensor_eigen.h"

using MessageLogger::Logger;
using MessageLogger::msg_log;


std::pair<MatrixOperator, RealDiagonalOperator>
SubspaceExpandBasis1(StateComponent& C, OperatorComponent const& H, StateComponent const& RightHam,
                     MixInfo const& Mix, KeepListType& KeepList,
		     std::set<QuantumNumbers::QuantumNumber> const& AddedQN,
		     StatesInfo const& States, TruncationInfo& Info,
                     StateComponent const& LeftHam)
{
   //TRACE(C);
   // truncate - FIXME: this is the s3e step
#if defined(SSC)
   MatrixOperator Lambda;
   SimpleStateComponent CX;
   std::tie(Lambda, CX) = ExpandBasis1_(C);
#else
   MatrixOperator Lambda = ExpandBasis1(C);
#endif

   //TRACE(Lambda);
   //TRACE(C);

   MatrixOperator Rho = scalar_prod(herm(Lambda), Lambda);
   //TRACE(Rho);
   if (Mix.MixFactor > 0)
   {
#if defined(SSC)
      StateComponent RH = contract_from_right(herm(H), CX, RightHam, herm(CX));
#else
      StateComponent RH = contract_from_right(herm(H), C, RightHam, herm(C));
#endif
      MatrixOperator RhoMix;
      MatrixOperator RhoL = scalar_prod(Lambda, herm(Lambda));

      // Skip the identity and the Hamiltonian
      for (unsigned i = 1; i < RH.size()-1; ++i)
      {
         double Prefactor = trace(triple_prod(herm(LeftHam[i]), RhoL, LeftHam[i])).real();
         if (Prefactor == 0)
            Prefactor = 1;
         //TRACE(Prefactor*triple_prod(herm(RH[i]), Rho, RH[i]));
         RhoMix += Prefactor * triple_prod(herm(RH[i]), Rho, RH[i]);
	 //TRACE(i)(Prefactor);
      }
      // check for a zero mixing term - can happen if there are no interactions that span
      // the current bond
      double RhoTrace = trace(RhoMix).real();
      if (RhoTrace != 0)
         Rho += (Mix.MixFactor / RhoTrace) * RhoMix;
   }
   if (Mix.RandomMixFactor > 0)
   {
      MatrixOperator RhoMix = MakeRandomMatrixOperator(Rho.Basis1(), Rho.Basis2());
      RhoMix = herm(RhoMix) * RhoMix;
      Rho += (Mix.RandomMixFactor / trace(RhoMix)) * RhoMix;
   }

   //   TRACE(Rho);

   DensityMatrix<MatrixOperator> DM(Rho);
   DM.DensityMatrixReport(std::cerr);
   DensityMatrix<MatrixOperator>::const_iterator DMPivot =
      TruncateFixTruncationErrorRelative(DM.begin(), DM.end(),
                                         States,
                                         Info);

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

   MatrixOperator UKeep = DM.ConstructTruncator(KeptStates.begin(), KeptStates.end());
   TRACE(UKeep);
   Lambda = Lambda * herm(UKeep);

   //TRACE(Lambda);

#if defined(SSC)
   C = UKeep*CX; //prod(U, CX);
#else
   C = prod(UKeep, C);
#endif

   MatrixOperator U, Vh;
   RealDiagonalOperator D;
   SVD_FullCols(Lambda, U, D, Vh);

   //TRACE(U)(D)(Vh);

   TRACE("QAZ")(U)(D)(Vh);

   C.debug_check_structure();

   C = prod(Vh, C);

   TRACE("WSX")(C);

   C.debug_check_structure();

   return std::make_pair(std::move(U), std::move(D));
}

// Apply subspace expansion / truncation on the right (C.Basis2()).
// Returns Lambda matrix (diagonal) and a unitary matrix
// Postcondition: C' Lambda' U' = C (up to truncation!)
std::pair<RealDiagonalOperator, MatrixOperator>
SubspaceExpandBasis2(StateComponent& C, OperatorComponent const& H, StateComponent const& LeftHam,
                     MixInfo const& Mix, KeepListType& KeepList,
		     std::set<QuantumNumbers::QuantumNumber> const& AddedQN,
		     StatesInfo const& States, TruncationInfo& Info,
                     StateComponent const& RightHam)
{
   //TRACE(C);

   // truncate - FIXME: this is the s3e step
   MatrixOperator Lambda = ExpandBasis2(C);

   MatrixOperator Rho = scalar_prod(Lambda, herm(Lambda));
   if (Mix.MixFactor > 0)
   {
      StateComponent LH = contract_from_left(H, herm(C), LeftHam, C);
      MatrixOperator RhoMix;

      MatrixOperator RhoR = scalar_prod(herm(Lambda), Lambda);

      for (unsigned i = 1; i < LH.size()-1; ++i)
      {
         double Prefactor = trace(triple_prod(herm(RightHam[i]), RhoR, RightHam[i])).real();
         if (Prefactor == 0)
            Prefactor = 1;
         RhoMix += Prefactor * triple_prod(LH[i], Rho, herm(LH[i]));
	 //	 TRACE(i)(Prefactor);
      }
      double RhoTrace = trace(RhoMix).real();
      if (RhoTrace != 0)
         Rho += (Mix.MixFactor / RhoTrace) * RhoMix;
   }
   if (Mix.RandomMixFactor > 0)
   {
      MatrixOperator RhoMix = MakeRandomMatrixOperator(Rho.Basis1(), Rho.Basis2());
      RhoMix = herm(RhoMix) * RhoMix;
      Rho += (Mix.RandomMixFactor / trace(RhoMix)) * RhoMix;
   }
   DensityMatrix<MatrixOperator> DM(Rho);
   DM.DensityMatrixReport(std::cerr);
   DensityMatrix<MatrixOperator>::const_iterator DMPivot =
      TruncateFixTruncationErrorRelative(DM.begin(), DM.end(),
                                         States,
                                         Info);
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

   MatrixOperator UKeep = DM.ConstructTruncator(KeptStates.begin(), KeptStates.end());

   Lambda = UKeep * Lambda;
   C = prod(C, herm(UKeep));

   MatrixOperator U, Vh;
   RealDiagonalOperator D;
   SVD_FullRows(Lambda, U, D, Vh);

   C = prod(C, U);

   return std::make_pair(std::move(D), std::move(Vh));
}



PStream::opstream& operator<<(PStream::opstream& out, DMRG const& d)
{
   return out << d.HamMatrices
              << d.Psi
	      << d.Site
              << d.Hamiltonian
	      << d.LeftStop
	      << d.RightStop

              << d.LastOverlap
              << d.IsPsiConverged
              << d.IsConvergedValid
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
              << d.SweepEntropy
              << d.SweepStartTime
              << d.SweepTruncatedEnergy
              << d.SweepEnergyError
              << d.SweepLastMixFactor

              << d.IterationNumMultiplies
              << d.IterationNumStates
              << d.IterationEnergy
              << d.IterationTruncation
              << d.IterationEntropy
              << d.IterationEnergyBeforeTrunc
              << d.IterationEnergyVec;
}

PStream::ipstream& operator>>(PStream::ipstream& in, DMRG& d)
{
   return in >> d.HamMatrices
             >> d.Psi
	     >> d.Site
             >> d.Hamiltonian
	     >> d.LeftStop
	     >> d.RightStop

             >> d.LastOverlap
             >> d.IsPsiConverged
             >> d.IsConvergedValid
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
             >> d.SweepEntropy
             >> d.SweepStartTime
             >> d.SweepTruncatedEnergy
             >> d.SweepEnergyError
             >> d.SweepLastMixFactor

             >> d.IterationNumMultiplies
             >> d.IterationNumStates
             >> d.IterationEnergy
             >> d.IterationTruncation
             >> d.IterationEntropy
             >> d.IterationEnergyBeforeTrunc
             >> d.IterationEnergyVec;

   // advance the C iterator
   int n = d.Site;
   d.C = d.Psi.begin();
   d.H = d.Hamiltonian.begin();
   while (n > 0)
   {
      ++d.C;
      ++d.H;
      --n;
   }
}

DMRG::DMRG(FiniteWavefunctionLeft const& Psi_, BasicTriangularMPO const& Ham_, int Verbose_)
   : Hamiltonian(Ham_),
     IsPsiConverged(false),
     NormalizeWavefunction(false),
     IsConvergedValid(false),
     MixUseEnvironment(false), UseDGKS(false), Solver_(LocalEigensolver::Solver::Lanczos),
     Verbose(Verbose_)

{
   // construct the HamMatrix elements
   if (Verbose > 0)
      std::cout << "Constructing Hamiltonian block operators..." << std::endl;
   H = Hamiltonian.begin();
   HamMatrices.push_left(Initial_E(Ham_, Psi_.Basis1()));
   TRACE(HamMatrices.left());
   for (FiniteWavefunctionLeft::const_mps_iterator I = Psi_.begin(); I != Psi_.end(); ++I)
   {
      if (Verbose > 1)
         std::cout << "site " << (HamMatrices.size_left()) << std::endl;
      HamMatrices.push_left(contract_from_left(*H, herm(*I), HamMatrices.left(), *I));
      Psi.push_back(*I);
      ++H;
   }
   HamMatrices.push_right(Initial_F(Ham_, Psi_.Basis2()));
   // Probably no need to incorporate lambda_r(), this is just a normalization
   // TODO: StateComponent * RealDiagonalOperator is not yet implemented
   Psi.set_back(prod(Psi.get_back().get(), MatrixOperator(Psi_.lambda_r().get())));

   // initialize to the right-most site
   HamMatrices.pop_left();
   C = Psi.end();
   --C;
   --H;
   Site = Psi.size()-1;

   LeftStop = 0;
   RightStop = Psi.size() - 1;

   // clear the global statistics
   TotalSweepNumber = 0;
   TotalSweepRecNumber = 0;
   TotalNumIterations = 0;
   TotalNumMultiplies = 0;

  this->debug_check_structure();
}

FiniteWavefunctionLeft
DMRG::Wavefunction() const
{
   return FiniteWavefunctionLeft::Construct(Psi);
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
      (*EnergyLog) << "#TotalSweepNumber #SweepNumber #Site "
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
   NormalizeWavefunction = Conf.Get("Normalize", true);
   Solver_.SetSolver(Conf.Get("Solver", "Lanczos"));
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
   NormalizeWavefunction = Conf.Get("Normalize", true);
   Solver_.SetSolver(Conf.Get("Solver", "Lanczos"));
   UseDGKS = Conf.Get("UseDGKS", false);

   // debugging; set the precision for cerr
   std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

   std::cout << "Convergence::OverlapTruncationRatio = " << ConvergenceOverlapTruncationScale << '\n';
   std::cout << "Convergence::OverlapDerivativeRatio = " << ConvergenceOverlapDifferenceOverlapScale << '\n';
   std::cout << "Convergence::TruncationCutoff = " << ConvergenceSweepTruncMin << '\n';
   std::cout << "Convergence::OverlapCutoff = " << ConvergenceOverlapMin << '\n';
   std::cout << "MixUseEnvironment = " << MixUseEnvironment << '\n';
   std::cout << "Normalize = " << NormalizeWavefunction << '\n';
   std::cout << "Solver = " << Solver_.GetSolverStr() << '\n';
   std::cout << "UseDGKS = " << UseDGKS << '\n';
}

void DMRG::StartIteration()
{
   IterationNumMultiplies = 0;
   IterationNumStates = 0;
   IterationEnergy = 0;
   IterationTruncation = 0;
   IterationEntropy = 0;
}

void DMRG::EndIteration()
{
   // recalculate the energy, as it will have changed after truncation
   if (IterationTruncation != 0)
      IterationEnergy = this->Energy();
   msg_log(1, "EnergyLog") << TotalSweepNumber << ' '
                           << TotalSweepRecNumber << ' '
                           << Site << ' '
                           << IterationNumStates << ' '
                           << IterationNumMultiplies << ' '
                           << IterationTruncation << ' '
                           << IterationEntropy << ' '
                           << IterationEnergy << '\n';

   ++SweepNumIterations;
   SweepNumMultiplies += IterationNumMultiplies;
   SweepSumStates += IterationNumStates;
   SweepMaxStates = std::max(SweepMaxStates, IterationNumStates);
   SweepEnergy = std::min(SweepEnergy, IterationEnergy.real());
   SweepTruncatedEnergy += IterationEnergyBeforeTrunc.real() - IterationEnergy.real();
   SweepTruncation += IterationTruncation;
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
   SweepEntropy = 0;
   SweepStartTime = ProcControl::GetCumulativeElapsedTime();
   SweepTruncatedEnergy = 0;
   IterationEnergyVec.clear();

   PsiPrevC = copy(*C);

   // Initialize the keep list
   KeepList.clear();

   IsConvergedValid = false;
}

double DMRG::FidelityLoss() const
{
#if 0
   return 2.0 * (1.0 - norm_frob(inner_prod(Psi.Center(), PsiPrevC))
                 / norm_frob(Psi.Center()));
#endif
   return 0;
}

void DMRG::EndSweep()
{
   TotalNumIterations += SweepNumIterations;
   TotalNumMultiplies += SweepNumMultiplies;

   SweepEnergyError = std::sqrt(statistics::variance(IterationEnergyVec.begin(),
                                                     IterationEnergyVec.end())
                                / SweepNumIterations);

   double Overlap = 2.0 * (1.0 - norm_frob(inner_prod(*C, PsiPrevC))
                           / norm_frob(*C));
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
                          << Overlap << ' '
                          << OverlapDifference << ' '
                          << SweepEnergy << ' '
                          << SweepLastMixFactor << ' '
                          << Converged
                          << std::endl;
   msg_log(1, "EnergyLog") << std::flush;
}

std::complex<double>
DMRG::Energy() const
{
   return inner_prod(HamMatrices.right(),
		     contract_from_left(*H, herm(HamMatrices.left()), *C, HamMatrices.right()));
}

std::complex<double>
DMRG::Solve()
{
   this->StartIteration();

#if 0
   // if the wavefunction is zero, add a random component
   double pNorm = norm_frob_sq(Psi.Center());
   if (pNorm < 1E-12)
   {
      Psi.Center() += MakeRandomMatrixOperator(Psi.Center().Basis1(), Psi.Center().Basis2());
      Psi.Center() *= 1.0 / norm_frob(Psi.Center());
   }
#endif

   Solver_.Solve(*C, HamMatrices.left(), *H, HamMatrices.right());

   DEBUG_TRACE("DMRG::Solve");
   StateComponent PsiOld = copy(*C);

   IterationEnergy = Solver_.LastEnergy();
   IterationNumMultiplies = Solver_.LastIter();

   IterationEnergyBeforeTrunc = IterationEnergy;

   return IterationEnergy;

   this->EndIteration();
}

bool DMRG::IsConverged() const
{
   return IsConvergedValid && IsPsiConverged;
}

TruncationInfo DMRG::TruncateAndShiftLeft(StatesInfo const& States)
{
   MatrixOperator U;
   RealDiagonalOperator Lambda;
   TruncationInfo Info;
   LinearWavefunction::const_iterator CNext = C;
   --CNext;

   TRACE(*C);

   std::tie(U, Lambda) = SubspaceExpandBasis1(*C, *H, HamMatrices.right(), MixingInfo,
					      KeepList, QuantumNumbersInBasis(CNext->LocalBasis()),
					      States, Info,
					      HamMatrices.left());

   if (Verbose > 1)
   {
      std::cerr << "Truncating left basis, states=" << Info.KeptStates() << '\n';
   }
   // update blocks

   TRACE(*C);

   HamMatrices.push_right(contract_from_right(herm(*H), *C, HamMatrices.right(), herm(*C)));

   TRACE(HamMatrices.right());

   HamMatrices.pop_left();

   // next site
   --Site;
   --H;
   --C;

   *C = prod(*C, U*Lambda);

   // normalize
   *C *= 1.0 / norm_frob(*C);

   IterationNumStates = Info.KeptStates();
   IterationTruncation += Info.TruncationError();
   IterationEntropy = std::max(IterationEntropy, Info.KeptEntropy());

   return Info;
}

TruncationInfo DMRG::TruncateAndShiftRight(StatesInfo const& States)
{
   // Truncate right
   RealDiagonalOperator Lambda;
   MatrixOperator U;
   TruncationInfo Info;
   LinearWavefunction::const_iterator CNext = C;
   ++CNext;
   std::tie(Lambda, U) = SubspaceExpandBasis2(*C, *H, HamMatrices.left(), MixingInfo,
					      KeepList, QuantumNumbersInBasis(CNext->LocalBasis()),
					      States, Info,
					      HamMatrices.right());
   if (Verbose > 1)
   {
      std::cerr << "Truncating right basis, states=" << Info.KeptStates() << '\n';
   }
   // update blocks
   HamMatrices.push_left(contract_from_left(*H, herm(*C), HamMatrices.left(), *C));
   HamMatrices.pop_right();

   // next site
   ++Site;
   ++H;
   ++C;

   *C = prod(Lambda*U, *C);

   // normalize
   *C *= 1.0 / norm_frob(*C);

   IterationNumStates = Info.KeptStates();
   IterationTruncation += Info.TruncationError();
   IterationEntropy = std::max(IterationEntropy, Info.KeptEntropy());

   return Info;
}

void DMRG::debug_check_structure() const
{
#if 0
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
#endif
}
