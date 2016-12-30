// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/krylov-time.h
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

#if !defined(KRYLOV_TIME_H_HHVERTHYEHYYOIUDLE89P)
#define KRYLOV_TIME_H_HHVERTHYEHYYOIUDLE89P

#include "simplekrylov.h"
#include "matrixproduct/mpstate.h"
#include "matrixproduct/mpoperator.h"
#include "matrixproduct/matrixproduct.h"
#include "matrixproduct/operatorstack.h"
#include "common/math_const.h"
#include <fstream>
#include <boost/shared_ptr.hpp>

struct TruncationInfo
{
   int m;
   double trunc;
   double entropy;
};

struct KrylovSolver
{
   typedef OperatorStack<MPStateComponent> SuperblockOperator;
   typedef OperatorStack<MatrixOperator>   TransformOperator;
   typedef MPWavefunction WavefunctionType;

   KrylovSolver() {}

   KrylovSolver(MPWavefunction const& Psi_,
                MPOperator const& H_,
                std::complex<double> Timestep_,
                int NumKrylovVectors_);

   void ShiftLeftAndExpand();
   void ShiftRightAndExpand();

   void ExpandLeft();
   void ExpandRight();

   int LeftSize() const { return Krylov.LeftSize(); }
   int RightSize() const { return Krylov.RightSize(); }

   // Updates the krylov basis.  This should be called at each iteration; it is the equivalent
   // of the Solve() function of other algorithms.
   void UpdateKrylovBasis();

   void MaximizeKrylovVectors() { Krylov.MaximizeKrylovVectors(); }
   void OrthogonalizeKrylovVectors() { Krylov.OrthogonalizeKrylovVectors(); }

   MPWavefunction& Wavefunction() { return Psi; }
   MPWavefunction const& Wavefunction() const { return Psi; }

   MPWavefunction const& LastKrylov() const { return Krylov.Krylov.back(); }

   void TruncateLeft(int MaxStates = DefaultMaxStates, double CFactor = 0);

   void TruncateRight(int MaxStates = DefaultMaxStates, double CFactor = 0);

   // Constructs Psi' = exp(Timestep * H) * Krylov[0]
   void AdvanceTimestep(int MaxStates);

   // debug sanity check that the wavefunction and operator basis agree
   void DebugCheckBasis() const;

   MPWavefunction Psi;            // initial wavefunction
   SimpleKrylov Krylov;            // the Krylov vectors

   QuantumNumber Ident;

   std::complex<double> Timestep;
};

PStream::opstream& operator<<(PStream::opstream& out, KrylovSolver const& s);
PStream::ipstream& operator>>(PStream::ipstream& in, KrylovSolver& s);

#endif
