// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp-algorithms/simplecg.h
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

#if !defined(DMRG_H_HHVERTHYEHYYOIUDLE89P)
#define DMRG_H_HHVERTHYEHYYOIUDLE89P

#include "matrixproduct/splitoperator.h"
#include "matrixproduct/centerwavefunction.h"
#include "matrixproduct/operatorstack.h"
#include "common/math_const.h"
#include "common/conflist.h"
#include <fstream>
#include <boost/shared_ptr.hpp>

double const Minus1OverPi = -math_const::r_1_pi;

struct Solver
{
   typedef OperatorStack<MPStateComponent> SuperblockOperator;
   typedef OperatorStack<MatrixOperator>   TransformOperator;
   typedef CenterWavefunction WavefunctionType;

   Solver() {}

   Solver(CenterWavefunction const& Psi_, SplitOperator const& Op_, CenterWavefunction const& Rhs_,
          double Freq_, double Broad_);

   virtual ~Solver();

   void StartSweep(bool IncrementSweepNumber = true,
                   double Broad = 0);  // prepare statistics for start of sweep

   void EndSweep();    // statistics for end of sweep

   void StartIteration();  // prepare statistics for start of iteration
   void EndIteration();    // statistics for end of iteration

   virtual void CreateLogFiles(std::string const& BasePath, ConfList const& Conf);
   virtual void RestoreLogFiles(std::string const& BasePath, ConfList const& Conf);

   // called by CreateLogFiles and RestoreLogFiles, loads the configuration options.
   virtual void ReadConfOptions(ConfList const& Conf);

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

   void ShiftLeft();
   void ShiftRight();

   void PrepareConvergenceTest() {}

   bool IsConverged() const { return false; }

   CenterWavefunction& Wavefunction();
   CenterWavefunction const& Wavefunction() const;

   CenterWavefunction WavefunctionAx() const;

   double ResidualNorm() const; // returns the residual norm (calculated, not stored)

   // returns the expectation value <x|A^2|x>
   virtual double ExpectationA2() const = 0;

   // returns the exact residual norm, given the supplied <x|A^2|x>
   virtual double ExactResidualNorm(double ExpectAt) const = 0;

   // returns the greens function via the DDMRG functional, given the supplied <x|A^2|x>
   virtual double Functional(double ExpectA2) const = 0;

   // returns the overlap < x | A | y >
   virtual std::complex<double> Overlap() const = 0;

   // returns the overlap <x|y>
   virtual std::complex<double> GreensFunction() const = 0;

   TruncationInfo TruncateLeft(StatesInfo const& SInfo, double CFactor);
   TruncationInfo TruncateRight(StatesInfo const& SInfo, double CFactor);

   // debug sanity check that the wavefunction and operator basis agree
   void DebugCheckBasis() const;

   CenterWavefunction x;             // the left hand side that we are solving for
   SplitOperator A;                 // the fixed operator
   CenterWavefunction y;             // the fixed right hand side
   MatrixOperator yprime;        // the right hand side in the mixed basis

   SuperblockOperator x_A_x;     // matrix elements   | x > A < x |
   TransformOperator x_y;        // matrix elements   | x > < y |

   QuantumNumber Ident;

   MatrixOperator LeftCorrection;  // Correction to the density matrix, left (basis 1)
   MatrixOperator RightCorrection;  // Correction to the density matrix, right (basis 2)

   double Frequency;
   double Broadening;

   // some global statistics
   int TotalSweepNumber;          // this is incremented every sweep
   int TotalSweepRecNumber;       // this is incremented conditionally
   int TotalNumIterations;
   int TotalNumMultiplies;

   // some statistics, for current sweep
   int SweepNumIterations;
   int SweepSumStates;            // sum of IterationNumStates
   int SweepMaxStates;
   int SweepNumMultiplies;        // number of mat-vec multiplies this sweep
   //   double SweepEnergy;            // lowest energy seen this sweep
   double SweepTruncation;        // cumulative truncation error this sweep
   double SweepEntropy;           // maximum entropy of the sweep
   double SweepGF_norm;
   double SweepGF_overlap;
   double SweepGF_real;

   // some statistics, for current iteration
   int IterationNumMultiplies;
   int IterationNumStates;
   double IterationTruncation;
   double IterationEntropy;
   double IterationGF_norm;
   double IterationGF_overlap;
   double IterationGF_real;
   double IterationSolverResid;

   // Log files are not serialized, but initialized by CreateLogFiles or
   // RestoreLogFiles
   boost::shared_ptr<std::ofstream> EnergyLog, DiagLog, SweepLog, CpuLog, DensityLog;

   // Following items are not persistent, but loaded from the configuration

   double LanczosMixFactor;   // the density matrix is rho_x + LanczosMixFactor * rho_y
   double AxMixFactor;        // Mix in a factor AxMixFactor * |Ax><Ax| (normalized)
   bool TruncateModifyScale;  // if true, use the modified weights in the density matrix
   bool FlushLogs;            // if true, flush logs after each output
};

PStream::opstream& operator<<(PStream::opstream& out, Solver const& s);
PStream::ipstream& operator>>(PStream::ipstream& in, Solver& s);

#endif
