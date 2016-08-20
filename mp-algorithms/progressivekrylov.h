// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/progressivekrylov.h
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

#if !defined(SIMPLEKRYLOV_HVKJH4378767T76RITFEW8REWH)
#define SIMPLEKRYLOV_HVKJH4378767T76RITFEW8REW
#include "matrixproduct/centerwavefunction.h"
#include "matrixproduct/splitoperator.h"
#include "matrixproduct/operatorstack.h"

using LinearAlgebra::Vector;
using LinearAlgebra::Matrix;

class KrylovSolver
{
   public:
      typedef OperatorStack<MPStateComponent> SuperblockOperator;
      typedef OperatorStack<MatrixOperator>   TransformOperator;
      typedef CenterWavefunction WavefunctionType;

      KrylovSolver();

      // The wavefunction must be normalized, prior to construction.
      KrylovSolver(SplitOperator const& H_, SplitOperator const& H2_, CenterWavefunction const& Psi);

      KrylovSolver(SplitOperator const& H_, SplitOperator const& H2_, CenterWavefunction const& Psi,
                   double Energy);

      KrylovSolver(SplitOperator const& H_, SplitOperator const& H2_, CenterWavefunction const& Psi,
                   double Energy, double ExpectationH2);

      void ShiftRightAndExpand(bool FullOrtho);
      void ShiftLeftAndExpand(bool FullOrtho);

      void ExpandLeft(bool FullOrtho);
      void ExpandRight(bool FullOrtho);

      int LeftSize() const { return Krylov.back().LeftSize(); }
      int RightSize() const { return Krylov.back().RightSize(); }

      // Adds another Krylov vector to the set, K is an initial guess vector
      void AddKrylovVector(CenterWavefunction const& K);

      // Once the next Krylov vector has converged, normalize it and calculate
      // matrix elements
      void FixKrylovVector(bool FullOrtho);

      Vector<std::complex<double> >
      Exponentiate(std::complex<double> Tau, double SubElement) const;

      Vector<std::complex<double> >
      Exponentiate(std::complex<double> Tau) const;

      void Solve(bool FullOrthogonalization);

      double Variance(bool FullOrtho) const;

   //double SubMatrixElement() const;

      void TruncateLeft(int MaxStates, double MixFactor, bool FullOrtho);
      void TruncateRight(int MaxStates, double MixFactor, bool FullOrtho);

      TruncationInfo TruncateLeft(int MinStates, int MaxStates, double MinTrunc, double MixFactor,
                                  bool FullMixing, bool FullOrtho);
      TruncationInfo TruncateRight(int MinStates, int MaxStates, double MinTrunc, double MixFactor,
                                   bool FullMixing, bool FullOrtho);

      // Constructs a linear combination of the Krylov vectors with coefficients given by n
      CenterWavefunction ConstructExpansion(LinearAlgebra::Vector<std::complex<double> > const& n,
                                            StatesInfo const& SInfo) const;

      // debug sanity check that the wavefunction and operator basis agree
      void DebugCheckBasis() const;

      std::vector<CenterWavefunction> Krylov;
      unsigned int NumConvergedKrylov;          // number of converged krylov vectors,
                                                // always either sizeof(Krylov) or sizeof(Krylov)-1
      SplitOperator H;
      SplitOperator H2;

      SuperblockOperator kn1_H_kn;              // matrix elements <k[n+1] | H | k[n] >

      Matrix<TransformOperator>       ki_kj;    // this is upper triangular, j > i

      Matrix<std::complex<double> >   sub_I;    // I in the Krylov basis
      Matrix<std::complex<double> >   sub_H;    // H in the Krylov basis, dimension is
                                                // number of Krylov vectors minus 1
      std::vector<double> sub_H2;               // Expectation values < k[i] | H^2 | k[i] >

      Matrix<std::complex<double> >   sub_L;    // Lower triangular matrix that orthogonalizes the Krylov vectors
      Matrix<std::complex<double> >   sub_Linv; // inverse of sub_L
      Matrix<std::complex<double> >   ortho_H;  // the Orthogonalized H

      QuantumNumber Ident;
};

PStream::opstream& operator<<(PStream::opstream& out, KrylovSolver const& s);
PStream::ipstream& operator>>(PStream::ipstream& in, KrylovSolver& s);

#endif
