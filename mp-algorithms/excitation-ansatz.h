// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/excitation-ansatz.h
//
// Copyright (C) 2022 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#if !defined(MPTOOLKIT_MP_ALGORITHMS_EXCITATION_ANSATZ_H)
#define MPTOOLKIT_MP_ALGORITHMS_EXCITATION_ANSATZ_H

#include "mpo/basic_triangular_mpo.h"
#include "mps/packunpack.h"
#include "wavefunction/ea.h"
#include "wavefunction/infinitewavefunctionleft.h"
#include "wavefunction/linearwavefunction.h"

// Struct to describe the effective Hamiltonian operator for the MPS excitation ansatz.
//
// HEff can be applied to a deque of MatrixOperators describing the X-matrices
// in the excitation ansatz, which are left-multiplied by the left null-space
// matrix to form the B-matrix such that
// PsiEA = \sum_n e^{ikn} ... AABAA...
// Each element of the deque describes the X-matrix at the position in the unit cell.
// At the moment, only the single-site EA is supported.
struct HEff
{
   HEff() {}

   // Initializer for the case where the A-matrices to the left and the right
   // of the excitation correspond to the same ground state Psi.
   HEff(InfiniteWavefunctionLeft const& Psi_, BasicTriangularMPO const& HamMPO_,
        QuantumNumbers::QuantumNumber const& Q_, ProductMPO const& StringOp_,
        double k, double GMRESTol_, int Verbose_);

   // Initializer for the case where PsiLeft and PsiRight are two DIFFERENT ground states.
   HEff(InfiniteWavefunctionLeft const& PsiLeft_, InfiniteWavefunctionLeft const& PsiRight_,
        BasicTriangularMPO const& HamMPO_, QuantumNumbers::QuantumNumber const& Q_,
        ProductMPO const& StringOp_, double k, double GMRESTol_, int Verbose_);

   // This function performs the part of the initialization common to both cases above.
   void Initialize();

   // Apply the effective Hamiltonian to XDeque.
   std::deque<MatrixOperator> operator()(std::deque<MatrixOperator> const& XDeque) const;

   // For a 2D system, if we specify the string operator describing Ty,
   // calculate the expectation value of Ty.
   std::complex<double> Ty(std::deque<MatrixOperator> const& XDeque) const;

   // Construct the B-matrices corresponding to the input X-matrices.
   std::deque<StateComponent> ConstructBDeque(std::deque<MatrixOperator> const& XDeque) const;

   // Construct the "triangular" MPS unit cell containing each B-matrix.
   LinearWavefunction ConstructPsiTri(std::deque<StateComponent> const& BDeque) const;

   // Generate a random initial state for a solver.
   // (This function is currently unused.)
   std::deque<MatrixOperator> InitialGuess() const;

   // Generate the objects to pack and unpack each X-matrix for the ARPACK wrapper.
   std::deque<PackMatrixOperator> PackInitialize() const;

   // Update the value of k.
   void SetK(double k);

   // Construct an EAWavefunction for the supplied X-matrices.
   EAWavefunction ConstructEAWavefunction(std::deque<MatrixOperator> XDeque) const;

   InfiniteWavefunctionLeft PsiLeft;
   InfiniteWavefunctionLeft PsiRight;
   BasicTriangularMPO HamMPO;
   QuantumNumbers::QuantumNumber Q;
   ProductMPO StringOp;
   double GMRESTol;
   int Verbose;
   std::complex<double> ExpIK;

   LinearWavefunction PsiLinearLeft, PsiLinearRight;
   StateComponent BlockHamL, BlockHamR;
   std::deque<StateComponent> BlockHamLDeque, BlockHamRDeque;
   std::deque<StateComponent> NullLeftDeque;
   // Eigenvectors of the left/right mixed transfer matrices.
   MatrixOperator RhoL, RhoR;
   // Eigenvectors of the left/right non-mixed (i.e. left/left and right/right)
   // transfer matrices with the Ty operator.
   MatrixOperator TyL, TyR;
   // Partially contracted versions of TyL and TyR.
   // Only needed when adding TyEff to HEff.
   //std::deque<StateComponent> TyLDeque, TyRDeque;
   // Left/right eigenvectors of the left/right mixed transfer matrices with
   // the Ty operator.
   MatrixOperator TyLRLeft, TyLRRight;
   MatrixOperator TyRLLeft, TyRLRight;
};

// A version of the HEff struct which can be used with the ARPACK wrapper.
struct PackHEff
{
   PackHEff(HEff H_);

   // Apply HEff to input vector.
   void operator()(std::complex<double> const* In_, std::complex<double>* Out_) const;

   // Convert raw array to deque of MatrixOperators.
   std::deque<MatrixOperator> unpack(std::complex<double> const* In_) const;

   // Convert deque of MatrixOperators to an array readable by the ARPACK wrapper.
   void pack(std::deque<MatrixOperator> XDeque, std::complex<double>* Out_) const;

   // Size of the raw arrays.
   int size() const
   {
      return Size;
   }

   std::deque<PackMatrixOperator> Pack;
   HEff H;
   int Size;
};

#endif
