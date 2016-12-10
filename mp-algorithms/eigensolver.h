// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/eigensolver.h
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

#if !defined(MPTOOLKIT_MP_ALGORITHMS_EIGENSOLVER_H)
#define MPTOOLKIT_MP_ALGORITHMS_EIGENSOLVER_H

#include "mps/state_component.h"
#include "mpo/operator_component.h"
#include "common/statistics.h"

class LocalEigensolver
{
   public:
      enum class Solver { InvalidSolver, Lanczos, Arnoldi, Davidson };

      static Solver SolverFromStr(std::string str);

      static std::string SolverStr(Solver s);

      LocalEigensolver();

      LocalEigensolver(Solver s);

      void SetInitialFidelity(int UnitCellSize, double f);

      // Apply the solver
      std::complex<double> Solve(StateComponent& C,
				 StateComponent const& LeftBlockHam,
				 OperatorComponent const& H,
				 StateComponent const& RightBlockHam);

      // Eigensolver parameters
      // Eigensolver tolerance is min(sqrt(AverageFidelity()) * FidelityScale, MaxTol)
      double FidelityScale;
      double MaxTol;
      double MinTol;

      // if EvolveDelta != 0 then do imaginary time evolution with this timestep instead of Lanczos
      double EvolveDelta;

      int MinIter; // Minimum number of iterations to perform (unless the eigensolver breaks down)
      int MaxIter; // Stop at this number, even if the eigensolver hasn't converged

      int Verbose;

      Solver GetSolver() const { return Solver_; }
      void SetSolver(Solver s);
      void SetSolver(std::string const& s) { this->SetSolver(SolverFromStr(s)); }

      // information on the state of the solver
      bool is_complex() const { return Solver_ == Solver::Arnoldi; }

      std::complex<double> LastEnergy() const { return LastEnergy_; }
      double LastEnergyReal() const { return LastEnergy_.real(); }
      double LastFidelity() const { return LastFidelity_; }
      double LastTol() const { return LastTol_; }
      double LastIter() const { return LastIter_; }
      double AverageFidelity() const { return FidelityAv_.value(); }

   private:
      Solver Solver_;
      // information on the last solver application
      double LastFidelity_;
      std::complex<double> LastEnergy_;
      double LastTol_;
      int LastIter_;

      // state information
      statistics::moving_exponential<double> FidelityAv_;
};

#endif
