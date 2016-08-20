// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/bicg.h
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

/*
  2005-10-12: NOTE: this header is now obsolete.
*/

#warning Obsolete

#if !defined(DMRG_H_HHVERTHYEHYYOIUDLE89P)
#define DMRG_H_HHVERTHYEHYYOIUDLE89P

#include "matrixproduct/mpstate.h"
#include "matrixproduct/mpoperator.h"
#include "matrixproduct/matrixproduct.h"
#include "operatorstack.h"

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

   MPWavefunction const& Wavefunction() const { return x; }

   int LeftSize() const { return x.LeftSize(); }
   int RightSize() const { return x.RightSize(); }

   // Calculates the wavefunction, using Iterations number
   // of Conjugate-Gradient steps
   double Solve(int Iterations);

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

   // returns the energy, calculated from a matrix-vector product
   double Energy();

   void Evaluate();

   TruncationInfo TruncateLeft(int MaxStates = DefaultMaxStates);

   TruncationInfo TruncateRight(int MaxStates = DefaultMaxStates);

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

};

#endif
