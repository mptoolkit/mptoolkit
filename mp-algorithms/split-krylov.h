// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/split-krylov.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#if !defined(SPLIT_KRYLOV_HVKJH4378767T76RITFEW8REWH)
#define SPLIT_KRYLOV_HVKJH4378767T76RITFEW8REW
#include "matrixproduct/mpstate.h"
#include "matrixproduct/mpoperator.h"
#include "matrixproduct/matrixproduct.h"
#include "matrixproduct/operatorstack.h"

using LinearAlgebra::Vector;
using LinearAlgebra::Matrix;

class SplitKrylov
{
   public:
      typedef OperatorStack<MPStateComponent> SuperblockOperator;
      typedef OperatorStack<MatrixOperator>   TransformOperator;
      typedef MPWavefunction WavefunctionType;

      SplitKrylov();

      SplitKrylov(MPOperator const& H_, std::vector<MPWavefunction> const& Krylov_);

      void ShiftRightAndExpand();
      void ShiftLeftAndExpand();

      void ExpandLeft();
      void ExpandRight();

      int LeftSize() const { return Krylov[0].LeftSize(); }
      int RightSize() const { return Krylov[0].RightSize(); }

      // Starting from Krylov[0], constructs an orthogonal Krylov basis
      // span{Krylov[i+1] = H * Krylov[i]}, setting the matrix elements
      // sub_H.
      void ConstructKrylovBasis();

      void ResetK0(MPWavefunction const& Psi);

      // Given a vector as an expansion of Krylov vectors |psi> = x[i] Krylov[i],
      // constructs the wavefunction |psi> = x[i] Krylov[i] in the Krylov[0] basis.
      MatrixOperator ConstuctWavefunctionFromRitz(Vector<std::complex<double> > const& x);

      // Given a vector as an expansion of Krylov vectors |psi> = x[i] Krylov[i],
      // constructs the reduced density matrix Tr_{L/R} |psi><psi| in the Krylov[0] basis.
      // The resulting density matrix is more accurate than the reduced density matrix
      // obtained from ConstructWavefunctionFromRitz, as fewer projections are done.
      MatrixOperator ConstructLeftDensityMatrixFromRitz(Vector<std::complex<double> > const& x);
      MatrixOperator ConstructRightDensityMatrixFromRitz(Vector<std::complex<double> > const& x);
 
      void TruncateLeft(double MixFactor, MatrixOperator Rho_k0L, 
			MatrixOperator Rho_k0R, int MaxStates);

      void TruncateRight(double MixFactor, MatrixOperator Rho_k0L, 
			 MatrixOperator Rho_k0R, int MaxStates);

      // The common part of TruncateLeft/Right
      void TruncateCommon(double MixFactor, MatrixOperator Rho_k0L, 
			  MatrixOperator Rho_k0R, int MaxStates, int Direction);

      // Returns the truncation operator U, given the density matrix Rho.
      // In principle, this does not depend on the left/right direction, or on i.  But
      // just in case we want to vary the truncation in some fancy way...
      MatrixOperator ConstructLeftTruncator(int i, MatrixOperator const& Rho, int MaxStates);
      MatrixOperator ConstructRightTruncator(int i, MatrixOperator const& Rho, int MaxStates);

      // Truncates the left basis of the i'th Krylov vector, using the operator U
      void TruncateKrylovLeft(std::size_t i, MatrixOperator const& U);

      // Truncates the right basis of the i'th Krylov vector, using the operator U
      void TruncateKrylovRight(std::size_t i, MatrixOperator const& U);

      // debug sanity check that the wavefunction and operator basis agree
      void DebugCheckBasis() const;

      std::vector<MPWavefunction> Krylov;
      MPOperator H;

      std::vector<SuperblockOperator> kn1_H_kn; // [0, NumKrylov()-1)
      Matrix<TransformOperator>       ki_kj;    // this is upper triangular, j > i

      Matrix<std::complex<double> >   sub_H;    // H in the Krylov basis, dimension is
                                                // number of Krylov vectors minus 1
      QuantumNumber Ident;
};

PStream::opstream& operator<<(PStream::opstream& out, SplitKrylov const& s);
PStream::ipstream& operator>>(PStream::ipstream& in, SplitKrylov& s);

#endif
