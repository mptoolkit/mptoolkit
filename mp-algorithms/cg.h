// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp-algorithms/cg.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

#if !defined(DMRG_H_HHVERTHYEHYYOIUDLE89P)
#define DMRG_H_HHVERTHYEHYYOIUDLE89P

#include "matrixproduct/mpstate.h"
#include "matrixproduct/mpoperator.h"
#include "matrixproduct/matrixproduct.h"
#include "matrixproduct/operatorstack.h"

struct TruncationInfo
{
   int m;
   double trunc;
   double entropy;
};

struct Solver
{
   typedef OperatorStack<MPStateComponent> SuperblockOperator;
   typedef OperatorStack<MatrixOperator>   TransformOperator;

   Solver(MPWavefunction const& Psi_, MPOperator const& Op_, MPWavefunction const& Rhs_);

   virtual ~Solver();

   MPWavefunction const& Wavefunction() const { return x; }

   int LeftSize() const { return x.LeftSize(); }
   int RightSize() const { return x.RightSize(); }

   // Calculates the wavefunction, using Iterations number
   // of Conjugate-Gradient steps.   Returns the norm of the residual.
   virtual double Solve(int Iterations) = 0;

   // Does a truncation and shifts the Center matrix to the right.
   // The new left basis is automatically expanded.
   void ShiftRightAndExpand();

   // Does a truncation and shifts the Center matrix to the left.
   // The new right basis is automatically expanded.
   void ShiftLeftAndExpand();

   // 'Expands' the left basis to cover all states.
   void ExpandLeft();

   // 'Expands' the right basis to cover all states.
   void ExpandRight();

   // returns the overlap < y | A | x >
   virtual std::complex<double> Overlap() const = 0;

   TruncationInfo TruncateLeft(int MaxStates = DefaultMaxStates, double CFactor = 0);

   TruncationInfo TruncateRight(int MaxStates = DefaultMaxStates, double CFactor = 0);

   // debug sanity check that the wavefunction and operator basis agree
   void DebugCheckBasis() const;

   MPWavefunction x;             // the left hand side that we are solving for
   MPOperator A;                 // the fixed operator
   MPWavefunction y;             // the fixed right hand side
   MPWavefunction yprime;        // the right hand side in the mixed basis

   SuperblockOperator yprime_A_x;  // matrix elements   | yprime > A < x |
   TransformOperator yprime_y;     // matrix elements   | yprime > < y |
   TransformOperator yprime_x;     // matrix elements   | yprime > < x |

   QuantumNumber Ident;

   MatrixOperator LeftCorrection;  // Correction to the density matrix, left (basis 1)
   MatrixOperator RightCorrection;  // Correction to the density matrix, right (basis 2)
};

#endif
