// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp-algorithms/solver-gmres.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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

#if !defined(SOLVER_GMRES_H_HDSFUIFHUIH89Y8943)
#define SOLVER_GMRES_H_HDSFUIFHUIH89Y8943

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/centerwavefunction.h"
#include "matrixproduct/splitoperator.h"
#include "pheap/pheap.h"
#include "simplecg.h"
#include "gmres.h"
#include "tensor/tensor_eigen.h"
#include "tensor/regularize.h"
#include <iostream>

class SolverGmres : public Solver
{
   public:
      SolverGmres() {}

      SolverGmres(CenterWavefunction const& Psi_, SplitOperator const& Op_,
                  SplitOperator const& Op2_, CenterWavefunction const& Rhs_,
                  double Freq_, double Broad_)
         : Solver(Psi_, Op_, Rhs_, Freq_, Broad_), ASquared(Op2_), Tol(0.0) {}

      virtual double Solve(int Iterations);

      virtual void ReadConfOptions(ConfList const& Conf);

      // returns the overlap <x|A|x>
      virtual std::complex<double> Overlap() const;

      virtual std::complex<double> GreensFunction() const;

      virtual double ExpectationA2() const;

      virtual double ExactResidualNorm(double ExpectA2) const;

      virtual double Functional(double ExpectA2) const;

      double NormSq() const { return norm_frob_sq(x.Center()); }

      SplitOperator ASquared;

      double Tol;  // tolerance of GMRES, read from configuration file
      int SubspaceSize; // GMRES krylov subspace size, read from config file
};

PStream::opstream& operator<<(PStream::opstream& out, SolverGmres const& s);
PStream::ipstream& operator>>(PStream::ipstream& in, SolverGmres& s);

#endif
