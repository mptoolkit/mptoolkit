// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/prodoptimizer.h
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

#if !defined(PRODOPTIMIZER_H_HHVERTHYEHYYOIUDLE89P)
#define PRODOPTIMIZER_H_HHVERTHYEHYYOIUDLE89P

#include "matrixproduct/splitoperator.h"
#include "matrixproduct/operatorstack.h"
#include "matrixproduct/centerwavefunction.h"

struct ProductOptimizer
{
   typedef OperatorStack<MPStateComponent> SuperblockOperator;

   ProductOptimizer(CenterWavefunction const& Psi_, SplitOperator const& A_, CenterWavefunction const& Rhs);

   ProductOptimizer(CenterWavefunction const& Psi_, std::vector<SplitOperator> const& A_,
                    std::vector<CenterWavefunction> const& Rhs_);

   CenterWavefunction const& Wavefunction() const { return Psi; }

   int LeftSize() const { return Psi.LeftSize(); }
   int RightSize() const { return Psi.RightSize(); }

   double Solve();

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

   TruncationInfo TruncateLeft(StatesInfo const& SInfo, double CFactor = 0);

   TruncationInfo TruncateRight(StatesInfo const& SInfo, double CFactor = 0);

   CenterWavefunction Psi;                 // we solve Psi = A * Rhs

   std::vector<SplitOperator> A;
   std::vector<CenterWavefunction> Rhs;

   std::vector<SuperblockOperator> Psi_A_Rhs;  // the matrix elements |Psi> A <Rhs|

   QuantumNumber Ident;

   void Initialize();

};

#endif
