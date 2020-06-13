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
#include "mp-algorithms/davidson.h"
#include <boost/algorithm/string.hpp>

#include "mps/packunpack.h"

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
   MPSMultiplyShift(StateComponent const& E_, OperatorComponent const& H_, StateComponent const& F_, std::complex<double> Energy_)
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
   std::complex<double> Energy;
};

LocalEigensolver::Solver
LocalEigensolver::SolverFromStr(std::string Str)
{
   Str = boost::to_lower_copy(Str);
   if (Str == "lanczos")
      return Solver::Lanczos;
   else if (Str == "arnoldi")
      return Solver::Arnoldi;
   else if (Str == "arnoldi-smallest")
      return Solver::ArnoldiSmallest;
   else if (Str == "arnoldi-lowest")
      return Solver::ArnoldiLowest;
   else if (Str == "davidson")
      return Solver::Davidson;
   else if (Str == "shift-invert")
      return Solver::ShiftInvert;
   else if (Str == "shift-invert-direct")
      return Solver::ShiftInvertDirect;
   else if (Str == "davidson")
      return Solver::Davidson;
   else if (Str == "davidson-target")
      return Solver::DavidsonTarget;
   else if (Str == "davidson-max-overlap")
      return Solver::DavidsonMaxOverlap;
   return Solver::InvalidSolver;
}

std::string
LocalEigensolver::SolverStr(LocalEigensolver::Solver s)
{
   if (s == Solver::Lanczos)
      return "lanczos";
   else if (s == Solver::Arnoldi)
      return "arnoldi";
   else if (s == Solver::Davidson)
      return "davidson";
   else if (s == Solver::ShiftInvert)
      return "shift-Invert";
   else if (s == Solver::ShiftInvertDirect)
      return "shift-invert-direct";
   else if (s == Solver::Davidson)
      return "davidson";
   else if (s == Solver::DavidsonTarget)
      return "davidson-target";
   else if (s == Solver::DavidsonMaxOverlap)
      return "davidson-max-overlap";
   return "Invalid Solver";
}

std::vector<std::string>
LocalEigensolver::EnumerateSolvers()
{
   std::vector<std::string> Result;
   for (int s = int(Solver::Lanczos); s <= int(Solver::LastSolver); ++s)
      Result.push_back(SolverStr(Solver(s)));
   return Result;
}

LocalEigensolver::LocalEigensolver(Solver s)
   : FidelityScale(0.1), MaxTol(1e-4), MinTol(1-10), MinIter(2), MaxIter(20), Verbose(0),
     Solver_(s), ShiftInvertEnergy(0), SubspaceSize(30), UsePreconditioning(false)
{
}

LocalEigensolver::LocalEigensolver()
   : LocalEigensolver(Solver::Lanczos)
{
}

void
LocalEigensolver::SetSolver(Solver s)
{
   Solver_ = s;
}

void
LocalEigensolver::SetShiftInvertEnergy(std::complex<double> E)
{
   ShiftInvertEnergy = E;
}

void
LocalEigensolver::SetSubspaceSize(int k)
{
   SubspaceSize = k;
}

void
LocalEigensolver::SetPreconditioning(bool p)
{
   UsePreconditioning = p;
}

void
LocalEigensolver::SetInitialFidelity(int UnitCellSize, double f)
{
   FidelityAv_ = statistics::moving_exponential<double>(exp(log(0.25)/UnitCellSize));
   FidelityAv_.push(f);
}

template <typename Func>
void
LinearSolveDirect(StateComponent& x, Func F, StateComponent const& Rhs, int Verbose = 0)
{
   LinearAlgebra::Matrix<double> HMat = LinearAlgebra::real(ConstructSuperOperator(F, x));

   PackStateComponent Pack(x);
   if (Verbose > 0)
   {
      std::cerr << "Linear solver dimension " << Pack.size() << '\n';
   }

   LinearAlgebra::Vector<std::complex<double>> v(Pack.size());
   Pack.pack(Rhs, v.data());

   LinearAlgebra::Matrix<double> vv(size(v),1);
   vv(LinearAlgebra::all,0) = LinearAlgebra::real(v);
   LinearAlgebra::Vector<std::complex<double>> xx = LinearAlgebra::LinearSolve(HMat, vv)(LinearAlgebra::all, 0);

   x = Pack.unpack(xx.data());
}

template <typename Func, typename Prec>
void
LinearSolve(StateComponent& x, Func F, Prec P, StateComponent const& Rhs, int k, int& MaxIter, double& Tol, int Verbose = 0)
{
   int m = k;     // krylov subspace size
   int iter = 0;   // total number of iterations performed

   double normb = norm_frob(P(Rhs));

   int IterThisRound = m*500;
   if (IterThisRound > MaxIter)
      IterThisRound = MaxIter;
   double tol = Tol;
   int Ret = GmRes(x, F, normb, Rhs, m, IterThisRound, tol, P, Verbose);
   iter += IterThisRound;

   while (Ret != 0 && iter < MaxIter)
   {
      // Attempt to avoid stagnation by increasing the number of iterations
      m += 10; // avoid stagnation
      if (Verbose > 1)
      {
	 std::cerr << "Refinement step, increasing m to " << m << '\n';
      }

      //      TRACE("Refinement step")(iter);
      // iterative refinement step
      StateComponent R = Rhs- F(x);
      StateComponent xRefine = R;
      IterThisRound = m*500;
      tol = Tol;
      Ret = GmRes(xRefine, F, normb, R, m, IterThisRound, tol, P, Verbose);

      iter += IterThisRound;

      if (Verbose > 1)
      {
	 double Resid = norm_frob(F(xRefine) - R) / normb;
	 std::cerr << "Residual of refined solver = " << Resid << '\n';
      }

      x += xRefine;

      if (Verbose > 1)
      {
	 double Resid = norm_frob(F(x) - Rhs) / normb;
	 std::cerr << "Residual after refinement step = " << Resid << '\n';
      }
   }
   Tol = tol;
   if (Ret != 0)
      Tol = -Tol;
   MaxIter = iter;
}

struct InverseDiagonalPrecondition
{
   InverseDiagonalPrecondition(StateComponent const& Diag_, std::complex<double> Energy_)
      : Diag(Diag_), Energy(Energy_)
   {
      //      TRACE(Diag);
      for (unsigned i = 0; i < Diag.size(); ++i)
      {
	 for (StateComponent::operator_type::iterator I = iterate(Diag[i]); I;  ++I)
	 {
	    for (StateComponent::operator_type::inner_iterator J = iterate(I); J; ++J)
	    {
	       for (auto II = iterate(*J); II; ++II)
	       {
		  for (auto JJ = iterate(II); JJ; ++JJ)
		  {
		     *JJ = (1.0 / (*JJ - Energy));
		  }
	       }
	    }
	 }
      }
      //TRACE(Diag)(Energy);
   }

   StateComponent operator()(StateComponent const& x) const
   {
      StateComponent Result(x);
      for (unsigned i = 0; i < x.size(); ++i)
      {
	 for (StateComponent::operator_type::iterator I = iterate(Result[i]); I;  ++I)
	 {
	    for (StateComponent::operator_type::inner_iterator J = iterate(I); J; ++J)
	    {
	       StateComponent::operator_type::const_inner_iterator d = iterate_at(Diag[i].data(), J.index1(), J.index2());
	       if (d)
	       {
		  *J = transform(*J, *d, LinearAlgebra::Multiplication<std::complex<double>, std::complex<double>>());
	       }
	    }
	 }
      }
      return Result;
   }

   StateComponent Diag;
   std::complex<double> Energy;

};

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
                   << " x " << C.Basis2().total_dimension() << ", requested Tol=" << LastTol_ << '\n';
      }
      if (Solver_ == Solver::Lanczos)
      {
	 LastEnergy_ = Lanczos(C, MPSMultiply(LeftBlockHam, H, RightBlockHam),
			       LastIter_, LastTol_, MinIter, Verbose-1);
      }
      else if (Solver_ == Solver::Arnoldi || Solver_ == Solver::ArnoldiSmallest)
      {
	 LastEnergy_ = LinearSolvers::Arnoldi(C, MPSMultiply(LeftBlockHam, H, RightBlockHam),
					      LastIter_, LastTol_,
					      LinearSolvers::SmallestMagnitude,
					      true, Verbose-1);
      }
      else if (Solver_ == Solver::ArnoldiLowest)
      {
	 LastEnergy_ = LinearSolvers::Arnoldi(C, MPSMultiply(LeftBlockHam, H, RightBlockHam),
					      LastIter_, LastTol_,
					      LinearSolvers::SmallestAlgebraicReal,
					      true, Verbose-1);
      }
      else if (Solver_ == Solver::ShiftInvert)
      {
	 StateComponent RHS = C;
	 if (UsePreconditioning)
	 {
	    if (Verbose > 2)
	    {
	       std::cerr << "Using diagonal preconditioning\n";
	    }
	    StateComponent D = operator_prod_inner_diagonal(H, LeftBlockHam, herm(RightBlockHam));

#if !defined(NDEBUG)
	    // debug check the diagonal
	    LinearAlgebra::Matrix<std::complex<double>> HMat = ConstructSuperOperator(MPSMultiplyShift(LeftBlockHam, H, RightBlockHam, ShiftInvertEnergy), D);
	    PackStateComponent Pack(D);
   	    LinearAlgebra::Vector<std::complex<double>> v(Pack.size());
	    Pack.pack(D, v.data());
	    LinearAlgebra::Vector<double> Diag = std::real(v);
	    for (int i = 0; i < size(Diag); ++i)
	    {
	       if (norm_frob(Diag[i] - HMat(i,i)) > 1E-10)
	       {
		  PANIC(i)(Diag[i])(HMat(i,i));
	       }
	    }
#endif

	    //	    TRACE(D);
	    LinearSolve(C, MPSMultiplyShift(LeftBlockHam, H, RightBlockHam, ShiftInvertEnergy),
			InverseDiagonalPrecondition(D, ShiftInvertEnergy),
			RHS, SubspaceSize, LastIter_, LastTol_, Verbose-1);
	    //TRACE(C);
	 }
	 else
	 {
	    LinearSolve(C, MPSMultiplyShift(LeftBlockHam, H, RightBlockHam, ShiftInvertEnergy),
			LinearAlgebra::Identity<StateComponent>(),
			RHS, SubspaceSize, LastIter_, LastTol_, Verbose-1);
	 }

	 // normalization
	 double CNorm = norm_frob(C);
	 //TRACE(CNorm);
	 // Adjust tol by CNorm to give a more realistic estimate
	 //LastTol_ /= CNorm;
	 C *= 1.0 / CNorm;
	 LastEnergy_ = inner_prod(C, MPSMultiply(LeftBlockHam, H, RightBlockHam)(C));
      }
      else if (Solver_ == Solver::ShiftInvertDirect)
      {
	 StateComponent RHS = C;
	 //TRACE(ShiftInvertEnergy);
	 LinearSolveDirect(C, MPSMultiplyShift(LeftBlockHam, H, RightBlockHam, ShiftInvertEnergy), RHS, Verbose-1);
	 LastTol_ = norm_frob(MPSMultiplyShift(LeftBlockHam, H, RightBlockHam, ShiftInvertEnergy)(C) - RHS);
	 // normalization
	 double CNorm = norm_frob(C);
	 //TRACE(CNorm);
	 // Adjust tol by CNorm to give a more realistic estimate
	 //LastTol_ /= CNorm;
	 C *= 1.0 / CNorm;
	 LastEnergy_ = inner_prod(C, MPSMultiply(LeftBlockHam, H, RightBlockHam)(C));
      }
      else if (Solver_ == Solver::DavidsonTarget)
      {
         CHECK(ShiftInvertEnergy.imag() == 0)("Davidson solver requires real energy");
         StateComponent Diagonal = operator_prod_inner_diagonal(H, LeftBlockHam, herm(RightBlockHam));
         LastEnergy_ = LinearSolvers::Davidson(C, Diagonal, MPSMultiply(LeftBlockHam, H, RightBlockHam),
                                               LinearSolvers::DavidsonMode::Target, LastIter_,
                                               Verbose-1, ShiftInvertEnergy.real());
      }
      else if (Solver_ == Solver::DavidsonMaxOverlap)
      {
         StateComponent Diagonal = operator_prod_inner_diagonal(H, LeftBlockHam, herm(RightBlockHam));
         LastEnergy_ = LinearSolvers::Davidson(C, Diagonal, MPSMultiply(LeftBlockHam, H, RightBlockHam),
                                               LinearSolvers::DavidsonMode::MaxOverlap, LastIter_,
                                               Verbose-1, ShiftInvertEnergy.real());
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

   LastFidelityLoss_ = std::max(1.0 - norm_frob(inner_prod(ROld, C)), 0.0);
   FidelityAv_.push(LastFidelityLoss_);
   return LastEnergy_;
}
