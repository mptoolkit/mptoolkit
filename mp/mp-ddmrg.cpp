// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-ddmrg.cpp
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
#include "matrixproduct/mpwavefunction.h"
#include "matrixproduct/mpoperator.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp-algorithms/ddmrg_functional.h"
#include "mp-algorithms/ddmrg.h"
#include "mp/copyright.h"
#include "tensor/tensor_eigen.h"
#include "tensor/regularize.h"
#include "interface/operator-parser.h"
#include <iostream>

class SolverDdmrg : public Solver
{
   public:
      SolverDdmrg(CenterWavefunction const& Psi_, SplitOperator const& Op_, 
                  CenterWavefunction const& Rhs_, double Freq, double Broad)
         : Solver(Psi_, Op_, Rhs_, Freq, Broad), Freq_(Freq), Broad_(Broad) {}

      virtual double Solve(int Iterations);

      virtual std::complex<double> Overlap() const;

      double Functional() const;

      double Freq_, Broad_;
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

struct InnerProd
{
   typedef std::complex<double> result_type;
   typedef MatrixOperator const& first_argument_type;
   typedef MatrixOperator const& second_argument_type;

   InnerProd(MatrixOperator const& Left, MatrixOperator const& Right) 
      : Left_(Left), Right_(Right) {}

   result_type operator()(MatrixOperator const& kx, MatrixOperator const& ky) const
   {
      return inner_prod(kx, triple_prod(Left_, ky, herm(Right_)));
   }

   MatrixOperator Left_, Right_;
};

double SolverDdmrg::Functional() const
{
   double Res = inner_prod(x.Center(), 
                           operator_prod(conj(A.Center()), 
                                         x_A_x.Left(), 
                                         x.Center(), 
                                         herm(x_A_x.Right()))).real();

   Res += 2.0 * Broad_ * 
      inner_prod(x.Center(), triple_prod(x_y.Left(), y.Center(), herm(x_y.Right()))).real();

   return Res;
}

double SolverDdmrg::Solve(int MaxIterations)
{
   //   TRACE(Psi.Center())(norm_frob(Psi.Center()));
   //   DEBUG_TRACE(Psi.Center().Basis1())(Psi.Center().Basis2());

   //   double Tol = 1E-10;
   //   int m = 20;

   double f = FunctionalMinimize(x.Center(), 
                                 SuperblockMultiply(conj(A.Center()),
                                                    x_A_x.Left(), 
                                                    x_A_x.Right()),
                                 Broad_ * triple_prod(x_y.Left(), y.Center(), herm(x_y.Right())),
                                 MaxIterations,
                                 //                      InnerProd(x_y.Left(), x_y.Right()),
                                 LinearAlgebra::Identity<MatrixOperator>());

   return f;
}

std::complex<double> SolverDdmrg::Overlap() const
{
   return inner_prod(operator_prod(A.Center(), 
                                   x_A_x.Left(), 
                                   x.Center(), 
                                   herm(x_A_x.Right())), 
                     Broad_ * triple_prod(x_y.Left(), y.Center(), herm(x_y.Right())));
}

int main(int argc, char** argv)
{
   if (argc < 7 || argc > 10)
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "DDMRG Solver\n";
      std::cerr << "usage: mp-ddmrg <Hamiltonian> <freq> <energy> <broad> <psi> <rhs> "
                << "[num-iter] [<maxstates>] [<correction>]\n";
      return 1;
   }

   TRACE(argc);

   pvalue_ptr<MPWavefunction> Psi = pheap::OpenPersistent(argv[5], 655360);
   OperatorList OpList;
   MPOperator H;
   std::tie(OpList, H) = ParseLatticeAndOperator(argv[1]);
   double Frequency = boost::lexical_cast<double>(argv[2]);
   double Energy =  boost::lexical_cast<double>(argv[3]);
   double Broadening = boost::lexical_cast<double>(argv[4]);
   pvalue_ptr<MPWavefunction> Rhs = pheap::ImportHeap(argv[6]);
   int NumIter = 4;
   if (argc >= 8) NumIter = boost::lexical_cast<int>(argv[7]);
   int MaxStates = 100;
   if (argc >= 9) MaxStates = boost::lexical_cast<int>(argv[8]);
   double CFactor = 0;
   if (argc >= 10) CFactor = boost::lexical_cast<double>(argv[9]);

   int MinStates = 0;
   double Trunc = 0;


   MPOperator Part = (Frequency+Energy) * OpList["I"] - H;
   MPOperator A = Part*Part + Broadening*Broadening*OpList["I"];

   SolverDdmrg solver(CenterWavefunction(*Psi), A, 
                      CenterWavefunction(*Rhs), Frequency, Broadening);

   std::cout.precision(14);

   CenterWavefunction xOld = solver.Wavefunction();
   //   MPWavefunction AxOld = solver.WavefunctionAx();

   bool First = false;
   for (int Sweeps = 0; Sweeps < 3; ++Sweeps)
   {

   {
      solver.ExpandLeft();
      double E = First ? 0.0 : solver.Solve(NumIter);
      TruncationInfo States = solver.TruncateLeft(MinStates, MaxStates, Trunc, CFactor);
      std::cout << '(' << solver.LeftSize() << ',' << solver.RightSize()
		<< ") " << E << ' ' << States.KeptStates() << '\n';
   }

   // sweep right
   while (solver.RightSize() > 1)
   {
      solver.ShiftRightAndExpand();
      //solver.ExpandRight();
      double E = First ? 0.0 : solver.Solve(NumIter);
      TruncationInfo States = solver.TruncateLeft(MinStates, MaxStates, Trunc, CFactor);
      std::cout << '(' << solver.LeftSize() << ',' << solver.RightSize()
		<< ") " << E << ' ' << States.KeptStates() << '\n';
   }
   First = false;


   //   std::cout << solver.ResidualNorm() << ' ' << solver.ExactResidualNorm() << '\n';
   std::cout //<< (difference(solver.Wavefunction(), xOld) / norm_2(xOld)) << ' '
      //             << (difference(solver.WavefunctionAx(), AxOld) / norm_2(AxOld)) << ' '
             << (overlap(solver.Wavefunction(), xOld) / norm_2_sq(xOld)) << '\n';
   xOld = solver.Wavefunction();
   //   AxOld = solver.WavefunctionAx();

   {
      solver.ExpandRight();
      double E = solver.Solve(NumIter);
      TruncationInfo States = solver.TruncateRight(MinStates, MaxStates, Trunc, CFactor);
      std::cout << '(' << solver.LeftSize() << ',' << solver.RightSize()
		<< ") " << E << ' ' << States.KeptStates() << '\n';

      TRACE(solver.Functional());
   }

   // sweep left
   while (solver.LeftSize() > 1)
   {
      solver.ShiftLeftAndExpand();
      //solver.ExpandLeft();
      double E = solver.Solve(NumIter);
      TruncationInfo States = solver.TruncateRight(MinStates, MaxStates, Trunc, CFactor);

      std::cout << '(' << solver.LeftSize() << ',' << solver.RightSize()
		<< ") " << E << ' ' << States.KeptStates() << '\n';
   }

   //   std::cout << solver.ResidualNorm() << ' ' << solver.ExactResidualNorm() << '\n';
   std::cout //<< (difference(solver.Wavefunction(), xOld) / norm_2(xOld)) << ' '
      //             << (difference(solver.WavefunctionAx(), AxOld) / norm_2(AxOld)) << ' '
             << (overlap(solver.Wavefunction(), xOld) / norm_2_sq(xOld)) << '\n';
   xOld = solver.Wavefunction();
   //   AxOld = solver.WavefunctionAx();

   }

   *Psi.mutate() = solver.Wavefunction().AsLinearWavefunction();

   pheap::ShutdownPersistent(Psi);
}
