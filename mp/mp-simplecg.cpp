// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-simplecg.cpp
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

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp-algorithms/simplecg.h"
#include "mp-algorithms/conjugategradient.h"
#include "mp/copyright.h"
#include <iostream>

class SolverHermitian : public Solver
{
   public:
      SolverHermitian(MPWavefunction const& Psi_, MPOperator const& Op_, MPWavefunction const& Rhs_)
         : Solver(Psi_, Op_, Rhs_) {}

      virtual double Solve(int Iterations);
      virtual double SolveLeft(int Iterations, int MaxStates);
      virtual double SolveRight(int Iterations, int MaxStates);

      virtual std::complex<double> Overlap() const;
};

struct SuperblockMultiply
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator const& argument_type;

   SuperblockMultiply(SimpleOperator const& Op_,
		      MPStateComponent const& Left_,
		      MPStateComponent const& Right_)
      : Op(Op_), Left(Left_), Right(Right_) {}

   MatrixOperator operator()(MatrixOperator const& Psi) const
   {
      return operator_prod(Op, Left, Psi, herm(Right));
   }

   SimpleOperator Op;
   MPStateComponent Left, Right;
};

struct SuperblockMultiplyLeft
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator const& argument_type;

   SuperblockMultiplyLeft(SimpleOperator const& Op_,
                          MPStateComponent const& Left_,
                          MPStateComponent const& Right_,
                          MatrixOperator const& Trunc_)
      : Op(Op_), Left(Left_), Right(Right_), Trunc(Trunc_) {}

   MatrixOperator operator()(MatrixOperator const& Psi) const
   {
      return prod(Trunc, operator_prod(Op, Left, Psi, herm(Right)), Trunc.TransformsAs());
   }

   SimpleOperator Op;
   MPStateComponent Left, Right;
   MatrixOperator Trunc;
};

struct SuperblockMultiplyRight
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator const& argument_type;

   SuperblockMultiplyRight(SimpleOperator const& Op_,
                           MPStateComponent const& Left_,
                           MPStateComponent const& Right_,
                           MatrixOperator const& Trunc_)
      : Op(Op_), Left(Left_), Right(Right_), Trunc(Trunc_) {}

   MatrixOperator operator()(MatrixOperator const& Psi) const
   {
      return prod(operator_prod(Op, Left, Psi, herm(Right)), Trunc, Trunc.TransformsAs());
   }

   SimpleOperator Op;
   MPStateComponent Left, Right;
   MatrixOperator Trunc;
};

struct RhsInner
{
   typedef std::complex<double> result_type;
   typedef MatrixOperator const& first_argument_type;
   typedef MatrixOperator const& second_argument_type;

   result_type operator()(MatrixOperator const& l, MatrixOperator const& r) const
   {
      return inner_prod(l, r);
   }
};

double SolverHermitian::Solve(int MaxIterations)
{
   //   TRACE(Psi.Center())(norm_frob(Psi.Center()));
   //   DEBUG_TRACE(Psi.Center().Basis1())(Psi.Center().Basis2());

   int Iter = MaxIterations;
   double Tol = 1E-10;

   ConjugateGradient(x.Center(), 
                     SuperblockMultiply(conj(A.Center()),
                                        x_A_x.Left(), 
                                        x_A_x.Right()),
                     yprime,
                     Iter, Tol,
                     LinearAlgebra::Identity<MatrixOperator>(),
                     RhsInner(),
                     RhsInner());
   return Tol;
}

double SolverHermitian::SolveLeft(int MaxIterations, int MaxStates)
{
   //   TRACE(Psi.Center())(norm_frob(Psi.Center()));
   //   DEBUG_TRACE(Psi.Center().Basis1())(Psi.Center().Basis2());

   int Iter = MaxIterations;
   double Tol = 1E-10;

   MatrixOperator Ax = operator_prod(conj(A.Center()), 
                                     x_A_x.Left(), 
                                     x.Center(), 
                                     herm(x_A_x.Right()));

   MatrixOperator Rho_Ax = scalar_prod(Ax, herm(Ax));

   DensityMatrix<MatrixOperator> DM(Rho_Ax * (1.0 / trace(Rho_Ax)));
   DensityMatrix<MatrixOperator>::const_iterator E = DM.begin();
   int Count = 0;
   while (E != DM.end() && Count < MaxStates)
   {
      ++E; ++Count;
   }
   MatrixOperator U = DM.ConstructTruncator(DM.begin(), E);
   U = scalar_prod(herm(U), U);
   MatrixOperator Yp = U * Ax;

   ConjugateGradient(x.Center(), 
                     SuperblockMultiply(conj(A.Center()),
                                        x_A_x.Left(), 
                                        x_A_x.Right()),
                     Yp,
                     Iter, Tol,
                     LinearAlgebra::Identity<MatrixOperator>(),
                     RhsInner(),
                     RhsInner());
   return Tol;
}

double SolverHermitian::SolveRight(int MaxIterations, int MaxStates)
{
   //   TRACE(Psi.Center())(norm_frob(Psi.Center()));
   //   DEBUG_TRACE(Psi.Center().Basis1())(Psi.Center().Basis2());

   int Iter = MaxIterations;
   double Tol = 1E-10;

   MatrixOperator Ax = operator_prod(conj(A.Center()), 
                                     x_A_x.Left(), 
                                     x.Center(), 
                                     herm(x_A_x.Right()));

   MatrixOperator Rho_Ax = scalar_prod(herm(Ax), Ax);
   DensityMatrix<MatrixOperator> DM(Rho_Ax * (1.0 / trace(Rho_Ax)));
   DensityMatrix<MatrixOperator>::const_iterator E = DM.begin();
   int Count = 0;
   while (E != DM.end() && Count < MaxStates)
   {
      ++E; ++Count;
   }
   MatrixOperator U = DM.ConstructTruncator(DM.begin(), E);
   U = scalar_prod(herm(U), U);
   MatrixOperator Yp = Ax * U;

   ConjugateGradient(x.Center(), 
                     SuperblockMultiply(conj(A.Center()),
                                        x_A_x.Left(), 
                                        x_A_x.Right()),
                     Yp,
                     Iter, Tol,
                     LinearAlgebra::Identity<MatrixOperator>(),
                     RhsInner(),
                     RhsInner());
   return Tol;
}

std::complex<double> SolverHermitian::Overlap() const
{
   return inner_prod(operator_prod(A.Center(), 
                                   x_A_x.Left(), 
                                   x.Center(), 
                                   herm(x_A_x.Right())), 
                     yprime);
}


int main(int argc, char** argv)
{
   if (argc < 5 || argc > 8)
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "Solver for complex Hermitian operators: A|x> = |y>\n";
      std::cerr << "usage: mp-simplecg <lattice> <operator> <psi> <rhs> "
                << "[num-iter] [<maxstates>] [<correction>]\n";
      return 1;
   }

   TRACE(argc);

   pvalue_ptr<MPWavefunction> Psi = pheap::OpenPersistent(argv[3], 655360);
   pvalue_ptr<OperatorList> System = pheap::ImportHeap(argv[1]);
   std::string Operator = argv[2];
   pvalue_ptr<MPWavefunction> Rhs = pheap::ImportHeap(argv[4]);
   int NumIter = 4;
   if (argc >= 6) NumIter = boost::lexical_cast<int>(argv[5]);
   int MaxStates = DefaultMaxStates;
   if (argc >= 7) MaxStates = boost::lexical_cast<int>(argv[6]);
   double CFactor = 0;
   if (argc >= 8) CFactor = boost::lexical_cast<double>(argv[7]);

   SolverHermitian solver(*Psi, (*System)[Operator], *Rhs);

   std::cout.precision(14);

   bool First = true;
   for (int Sweeps = 0; Sweeps < 3; ++Sweeps)
   {

   {
      solver.ExpandLeft();
      TRACE(solver.ResidualNorm());
      double E = First ? 0.0 : solver.Solve(NumIter);
      TRACE(solver.ResidualNorm());
      TruncationInfo States = solver.TruncateLeft(MaxStates, CFactor);
      TRACE(solver.ResidualNorm());
      std::cout << '(' << solver.LeftSize() << ',' << solver.RightSize()
		<< ") " << E << ' ' << States.m << std::endl;
   }

   // sweep right
   while (solver.RightSize() > 1)
   {

      if (solver.RightSize() == 22 && !First)
      {
         
         for (int i = 0; i < 10; ++i)
         {
            TRACE(solver.ResidualNorm());
            solver.ShiftRightAndExpand();
            //            solver.ExpanfRight();
            TRACE(solver.ResidualNorm());
            solver.ExpandLeft();
            TRACE(solver.ResidualNorm());
            double E = First ? 0.0 : solver.Solve(NumIter);
            TRACE(solver.ResidualNorm());
            TruncationInfo States = solver.TruncateLeft(MaxStates, CFactor);
            TRACE(solver.ResidualNorm());
            
            std::cout << '(' << solver.LeftSize() << ',' << solver.RightSize()
                      << ") " << E << ' ' << States.m << std::endl;


#if 0
            solver.ExpandLeft();
            TRACE(solver.ResidualNorm());
            E = solver.Solve(NumIter);
            TRACE(solver.ResidualNorm());
            std::cout << '(' << solver.LeftSize() << ',' << solver.RightSize()
                      << ") " << E << ' ' << States.m << ' ' << States.trunc << std::endl;

            solver.ExpandLeft();
            TRACE(solver.ResidualNorm());
            E = solver.Solve(NumIter);
            TRACE(solver.ResidualNorm());
            States = solver.TruncateLeft(MaxStates, CFactor);
            TRACE(solver.ResidualNorm());
            States = solver.TruncateLeft(MaxStates, CFactor);
            TRACE(solver.ResidualNorm());
            solver.ExpandLeft();
            TRACE(solver.ResidualNorm());
            States = solver.TruncateLeft(MaxStates, CFactor);
            TRACE(solver.ResidualNorm());
            std::cout << '(' << solver.LeftSize() << ',' << solver.RightSize()
                      << ") " << E << ' ' << States.m << ' ' << States.trunc << std::endl;

#endif
            solver.ShiftLeftAndExpand();
            //      solver.ExpandLeft();
            TRACE(solver.ResidualNorm());
            E = solver.Solve(NumIter);
            TRACE(solver.ResidualNorm());
            States = solver.TruncateRight(MaxStates, CFactor);
            TRACE(solver.ResidualNorm());

            std::cout << '(' << solver.LeftSize() << ',' << solver.RightSize()
                      << ") " << E << ' ' << States.m << ' ' << States.trunc << std::endl;
         }
      }

      TRACE(solver.ResidualNorm());
      solver.ShiftRight();
      TRACE(solver.ResidualNorm());
      solver.ShiftLeft();
      TRACE(solver.ResidualNorm());
      solver.ShiftRight();
      TRACE(solver.ResidualNorm());
      solver.ExpandLeft();
      //      solver.ExpandRight();
      TRACE(solver.ResidualNorm());
      double E = First ? 0.0 : solver.Solve(NumIter);
      TRACE(solver.ResidualNorm());
      E = First ? 0.0 : solver.SolveLeft(NumIter, MaxStates);
      TRACE(solver.ResidualNorm());

      if (solver.LeftSize() == 10)
      {
         TRACE("Truncating with max 10000 states:");
         TruncationInfo States = solver.TruncateLeft(10000);
         TRACE(solver.ResidualNorm());
      }

      TruncationInfo States = solver.TruncateLeft(MaxStates, CFactor);
      TRACE(solver.ResidualNorm());

      std::cout << '(' << solver.LeftSize() << ',' << solver.RightSize()
		<< ") " << E << ' ' << States.m << ' ' << States.trunc << std::endl;
   }
   First = false;

   {
      solver.ExpandRight();
      double E = solver.Solve(NumIter);
      TruncationInfo States = solver.TruncateRight(MaxStates, CFactor);
      std::cout << '(' << solver.LeftSize() << ',' << solver.RightSize()
		<< ") " << E << ' ' << States.m << ' ' << States.trunc << std::endl;
   }

   // sweep left
   while (solver.LeftSize() > 1)
   {
      solver.ShiftLeftAndExpand();
      //      solver.ExpandLeft();
      double E = solver.Solve(NumIter);
      TRACE(solver.ResidualNorm());
      E = solver.SolveRight(NumIter, MaxStates);
      TRACE(solver.ResidualNorm());
      TruncationInfo States = solver.TruncateRight(MaxStates, CFactor);

      std::cout << '(' << solver.LeftSize() << ',' << solver.RightSize()
		<< ") " << E << ' ' << States.m << ' ' << States.trunc << std::endl;
   }

   }

   *Psi.mutate() = solver.Wavefunction();

   pheap::ShutdownPersistent(Psi);
}
