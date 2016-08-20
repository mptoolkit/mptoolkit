// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-simplecg-sym.cpp
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

class SolverSymmetric : public Solver
{
   public:
      SolverSymmetric(MPWavefunction const& Psi_, MPOperator const& Op_, MPWavefunction const& Rhs_)
         : Solver(Psi_, Op_, Rhs_) {}

      virtual double Solve(int Iterations);

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

struct RhsInner
{
   typedef std::complex<double> result_type;
   typedef MatrixOperator const& first_argument_type;
   typedef MatrixOperator const& second_argument_type;

   result_type operator()(MatrixOperator const& l, MatrixOperator const& r) const
   {
      return conj(inner_prod(l, conj(r)));
   }
};

struct LhsInner
{
   typedef std::complex<double> result_type;
   typedef MatrixOperator const& first_argument_type;
   typedef MatrixOperator const& second_argument_type;

   result_type operator()(MatrixOperator const& l, MatrixOperator const& r) const
   {
      return conj(inner_prod(l, conj(r)));
   }
};

double SolverSymmetric::Solve(int MaxIterations)
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
                     LhsInner());

   return Tol;
}

std::complex<double> SolverSymmetric::Overlap() const
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
      std::cerr << "Solver for complex Symmetric operators: A|x> = |y>\n";
      std::cerr << "usage: mp-cg-sym <lattice> <operator> <psi> <rhs> "
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

   SolverSymmetric solver(*Psi, (*System)[Operator], *Rhs);

   std::cout.precision(14);

   bool First = false;
   for (int Sweeps = 0; Sweeps < 3; ++Sweeps)
   {

   {
      solver.ExpandLeft();
      double E = First ? 0.0 : solver.Solve(NumIter);
      TruncationInfo States = solver.TruncateLeft(MaxStates, CFactor);
      std::cout << '(' << solver.LeftSize() << ',' << solver.RightSize()
                << ") " << E << ' ' << States.m << '\n';
   }

   // sweep right
   while (solver.RightSize() > 1)
   {
      solver.ShiftRightAndExpand();
      solver.ExpandRight();
      double E = First ? 0.0 : solver.Solve(NumIter);
      TruncationInfo States = solver.TruncateLeft(MaxStates, CFactor);

      std::cout << '(' << solver.LeftSize() << ',' << solver.RightSize()
                << ") " << E << ' ' << States.m << '\n';
   }
   First = false;

   {
      solver.ExpandRight();
      double E = solver.Solve(NumIter);
      TruncationInfo States = solver.TruncateRight(MaxStates, CFactor);
      std::cout << '(' << solver.LeftSize() << ',' << solver.RightSize()
                << ") " << E << ' ' << States.m << '\n';
   }

   // sweep left
   while (solver.LeftSize() > 1)
   {
      solver.ShiftLeftAndExpand();
      solver.ExpandLeft();
      double E = solver.Solve(NumIter);
      TruncationInfo States = solver.TruncateRight(MaxStates, CFactor);

      std::cout << '(' << solver.LeftSize() << ',' << solver.RightSize()
                << ") " << E << ' ' << States.m << '\n';
   }

   }

   *Psi.mutate() = solver.Wavefunction();

   pheap::ShutdownPersistent(Psi);
}
