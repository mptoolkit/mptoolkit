// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/krylovloop.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER
/* -*- C++ -*- $Id$
 */

#if !defined(KRYLOVLOOP_H_CDHJH43Y57Y78Y87H78YT875YORE)
#define KRYLOVLOOP_H_CDHJH43Y57Y78Y87H78YT875YORE

#include "mp-algorithms/progressivekrylov.h"

class KrylovLoop
{
   public:
      typedef KrylovSolver SolverType;
      typedef SolverType::WavefunctionType WavefunctionType;

      KrylovLoop() {}

      KrylovLoop(SplitOperator const& H_, SplitOperator const& H2_, WavefunctionType const& Psi) 
         : Solver_(H_, H2_, (1.0 / norm_frob(Psi.Center())) * Psi), PsiNorm_(norm_frob(Psi.Center())) {}

      KrylovLoop(SolverType const& Solver);

      double SweepRight(double MinTrunc, bool FullMixing);
      double SweepLeft(double MinTrunc, bool FullMixing);

      // Adds another Krylov vector, optionally using K as an initial guess
      void AddKrylovVector(WavefunctionType const& K);
      void AddKrylovVector();

      // Sweeps the just added Krylov vector until convergence.  Note the parameter
      // here is the Variance, not the error!
      void EvolveNextKrylov(double RequiredVariance);

      double GetTimestep(double GuessTimestep) const;

      // This does everything, starting from an empty Krylov basis,
      // building a fixed size krylov basis, using GuessTimestep as a guide
      // for the error terms.
      void ConstructKrylovBasis(int NumKrylov, double GuessTimestep);

      // This does everything, starting from an empty Krylov basis,
      // constructing a Krylov basis until the error is small enough.
      void ConstructKrylovBasis(std::complex<double> Timestep);

      WavefunctionType Wavefunction() const;

      static double TimestepFromHochbruckBound(int NumKrylov, double HNorm, double EBound);

      Matrix<std::complex<double> > SubspaceIdentity() const { return Solver_.sub_I; }
      Matrix<std::complex<double> > SubspaceH() const { return Solver_.sub_H; }

      void Evolve(double Timestep, StatesInfo const& SInfo);
      void Evolve(std::complex<double> Timestep, StatesInfo const& SInfo);

      void Evolve(double Timestep, StatesInfo const& SInfo, SplitOperator const& NewH);

   // returns the norm of the wavefunction
      double Norm() const { return PsiNorm_; }

   //private:
      SolverType Solver_;

      double PsiNorm_;
      double LastTimestep_;

      double ErrorBound_;  // bound on the error ||psi(exact) - psi||^2 per unit timestep
      // FIXME: this has a different interpretation depending on which version of ConstructKrylovBasis is used.
      std::complex<double> TimeDirection_;  // this has unit norm, it is the complex direction to evolve in

      int MinStates_;
      int MaxStates_;
      bool TwoSite_;
      bool FullOrtho_;
      double MixFactor_;
      double RealTime_;
      double Beta_;
      double ErrorScaleFactor_;
      int LastNumKrylov_;
      bool DoReductionSweep_;
      double KrylovCutoffFactor_;
      int MaxSweeps_;
      
      Vector<std::complex<double> > OldCoeffVec_;
};

#endif
