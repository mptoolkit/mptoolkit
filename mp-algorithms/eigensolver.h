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
      // FIXME: this could use common/namedenum.h
      enum class Solver { InvalidSolver, Lanczos, Arnoldi, ArnoldiSmallest, ArnoldiLowest, ShiftInvert, ShiftInvertDirect,
                             Davidson, DavidsonTarget, DavidsonMaxOverlap,
                             LastSolver = DavidsonMaxOverlap};

      static Solver SolverFromStr(std::string str);

      static std::string SolverStr(Solver s);

      static std::vector<std::string> EnumerateSolvers();

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
      std::string GetSolverStr() const { return SolverStr(Solver_); }
      void SetSolver(Solver s);
      void SetSolver(std::string const& s) { this->SetSolver(SolverFromStr(s)); }

      // For solver-specific parameters, we probably should provide a generic format.
      // Since we only have one solver-specific parameter we just have a special case for now.
      void SetShiftInvertEnergy(std::complex<double> E);

      // subspace size, used by shift-invert solver
      void SetSubspaceSize(int k);

      // set whether or not to use preconditioning in the shift-invert solver
      void SetPreconditioning(bool Pre);

      // A better approach is a function to set a solver parameter, eg
      // void SetSolverParameter(std::string const& s);

      // information on the state of the solver
      bool is_complex() const { return Solver_ == Solver::Arnoldi || Solver_ == Solver::ArnoldiLowest; }

      std::complex<double> LastEnergy() const { return LastEnergy_; }
      double LastEnergyReal() const { return LastEnergy_.real(); }
      double LastFidelityLoss() const { return LastFidelityLoss_; }
      double LastTol() const { return LastTol_; }
      double LastIter() const { return LastIter_; }
      double AverageFidelity() const { return FidelityAv_.value(); }

   private:
      Solver Solver_;
      // information on the last solver application
      double LastFidelityLoss_;
      std::complex<double> LastEnergy_;
      double LastTol_;
      int LastIter_;
      std::complex<double> ShiftInvertEnergy;
      int SubspaceSize;
      bool UsePreconditioning;

      // state information
      statistics::moving_exponential<double> FidelityAv_;
};

#endif
