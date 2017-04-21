// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/eigensolver.cpp
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

#include "eigensolver.h"
#include "mp-algorithms/lanczos.h"
#include "mp-algorithms/arnoldi.h"
#include "mp-algorithms/gmres.h"
#include <boost/algorithm/string.hpp>

struct MPSMultiply
{
   MPSMultiply(StateComponent const& E_, OperatorComponent const& H_, StateComponent const& F_)
      : E(E_), H(H_), F(F_)
   {
   }

   StateComponent operator()(StateComponent const& Psi) const
   {
      StateComponent R = operator_prod_inner(H, E, Psi, herm(F));
      return R;
   }

   StateComponent const& E;
   OperatorComponent const& H;
   StateComponent const& F;
};

struct MPSMultiplyShift
{
   MPSMultiplyShift(StateComponent const& E_, OperatorComponent const& H_, StateComponent const& F_, double Energy_)
      : E(E_), H(H_), F(F_), Energy(Energy_)
   {
   }

   StateComponent operator()(StateComponent const& Psi) const
   {
      StateComponent R = operator_prod_inner(H, E, Psi, herm(F)) - Energy*Psi;
      return R;
   }

   StateComponent const& E;
   OperatorComponent const& H;
   StateComponent const& F;
   double Energy;
};

LocalEigensolver::Solver
LocalEigensolver::SolverFromStr(std::string Str)
{
   Str = boost::to_lower_copy(Str);
   if (Str == "lanczos")
      return Solver::Lanczos;
   else if (Str == "arnoldi")
      return Solver::Arnoldi;
   else if (Str == "davidson")
      return Solver::Davidson;
   else if (Str == "shift-invert")
      return Solver::ShiftInvert;
   return Solver::InvalidSolver;
}

std::string
LocalEigensolver::SolverStr(LocalEigensolver::Solver s)
{
   if (s == Solver::Lanczos)
      return "Lanczos";
   else if (s == Solver::Arnoldi)
      return "Arnoldi";
   else if (s == Solver::Davidson)
      return "Davidson";
   else if (s == Solver::ShiftInvert)
      return "Shift-Invert";
   return "Invalid Solver";
}

LocalEigensolver::LocalEigensolver()
   : FidelityScale(0.1), MaxTol(1e-4), MinTol(1-10), MinIter(2), MaxIter(20), Verbose(0),
     Solver_(Solver::Lanczos), ShiftInvertEnergy(0)
{
}

LocalEigensolver::LocalEigensolver(Solver s)
   : FidelityScale(0.1), MaxTol(1e-4), MinTol(1-10), MinIter(2), MaxIter(20), Verbose(0),
     Solver_(s), ShiftInvertEnergy(0)
{
}

void
LocalEigensolver::SetSolver(Solver s)
{
   Solver_ = s;
}

void
LocalEigensolver::SetShiftInvertEnergy(double E)
{
   ShiftInvertEnergy = E;
}

void
LocalEigensolver::SetInitialFidelity(int UnitCellSize, double f)
{
   FidelityAv_ = statistics::moving_exponential<double>(exp(log(0.25)/UnitCellSize));
   FidelityAv_.push(f);
}

double CNorm = 1;
double ActualEnergy = 1;
int SubspaceSize = 30;

std::complex<double>
LocalEigensolver::Solve(StateComponent& C,
                        StateComponent const& LeftBlockHam,
                        OperatorComponent const& H,
                        StateComponent const& RightBlockHam)
{
   DEBUG_CHECK_EQUAL(C.Basis1(), LeftBlockHam.Basis2());
   DEBUG_CHECK_EQUAL(C.Basis2(), RightBlockHam.Basis1());
   DEBUG_CHECK_EQUAL(LeftBlockHam.LocalBasis(), H.Basis1());
   DEBUG_CHECK_EQUAL(RightBlockHam.LocalBasis(), H.Basis2());
   DEBUG_CHECK_EQUAL(C.LocalBasis(), H.LocalBasis2());

   StateComponent ROld = C;

   if (EvolveDelta == 0.0)
   {
      LastTol_ = std::min(std::sqrt(this->AverageFidelity()) * FidelityScale, MaxTol);
      LastTol_ = std::max(LastTol_, MinTol);
      //LastTol_ = std::min(this->AverageFidelity() * FidelityScale, MaxTol);
      LastIter_ = MaxIter;
      if (Verbose > 2)
      {
         std::cerr << "Starting eigensolver.  Initial guess vector has dimensions "
                   << C.Basis1().total_dimension() << " x " << C.LocalBasis().size()
                   << " x " << C.Basis2().total_dimension() << '\n';
      }
      if (Solver_ == Solver::Lanczos)
      {
	 LastEnergy_ = Lanczos(C, MPSMultiply(LeftBlockHam, H, RightBlockHam),
			       LastIter_, LastTol_, MinIter, Verbose-1);
      }
      else if (Solver_ == Solver::Arnoldi)
      {
	 LastEnergy_ = LinearSolvers::Arnoldi(C, MPSMultiply(LeftBlockHam, H, RightBlockHam),
					      LastIter_, LastTol_, 
					      LinearSolvers::SmallestMagnitude,
					      true, Verbose-1);
      }
      else if (Solver_ == Solver::ShiftInvert)
      {
	 StateComponent RHS = C;
	 // scale the initial state using the previous norm
	 //C *= CNorm;
	 //TRACE(norm_frob(C));
	 //ActualEnergy = inner_prod(C, MPSMultiply(LeftBlockHam, H, RightBlockHam)(C)).real();
	 //C *= 1.0 / (ActualEnergy - ShiftInvertEnergy);
	 C *= 0;
         GmRes(C, MPSMultiplyShift(LeftBlockHam, H, RightBlockHam, ShiftInvertEnergy),
	       RHS, SubspaceSize, LastIter_, LastTol_,
	       LinearAlgebra::Identity<StateComponent>(), Verbose-1);
	 // normalization
	 CNorm = norm_frob(C);
	 TRACE(CNorm);
	 // Adjust tol by CNorm to give a more realistic estimate
	 //LastTol_ /= CNorm;
	 C *= 1.0 / CNorm;
	 ActualEnergy = inner_prod(C, MPSMultiply(LeftBlockHam, H, RightBlockHam)(C)).real();
	 //TRACE(ActualEnergy);
	 LastEnergy_ = ActualEnergy; //ShiftInvertEnergy;
      }
      else
      {
	 PANIC("Unsupported solver")(SolverStr(Solver_));
      }
   }
   else
   {
      C = operator_prod_inner(H, LeftBlockHam, ROld, herm(RightBlockHam));
      LastEnergy_ = inner_prod(ROld, C).real();
      C = ROld - EvolveDelta * C; // imaginary time evolution step
      C *= 1.0 / norm_frob(C);    // normalize
      LastIter_ = 1;
      LastTol_ = 0.0;
   }

   LastFidelity_ = std::max(1.0 - norm_frob(inner_prod(ROld, C)), 0.0);
   FidelityAv_.push(LastFidelity_);
   return LastEnergy_;
}
