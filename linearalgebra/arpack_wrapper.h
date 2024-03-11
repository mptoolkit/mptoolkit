// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/arpack_wrapper.h
//
// Copyright (C) 2007-2022 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2022 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#if !defined(ARPACK_WRAPPER_H_DSHUY7R894789LX)
#define ARPACK_WRAPPER_H_DSHUY7R894789LX

#include "common/trace.h"
#include "linearalgebra/vector.h"
#include <vector>
#include <complex>

namespace LinearAlgebra
{

// Use ARPACK to find eigenvalues and optionally eigenvectors
// of a complex non-hermitian matrix.  Mult is a functor that
// takes two arguments, an input vector of n complex numbers, and
// writes the result of the matrix-vector multiply to the second output argument.
// TODO: this could easily be generalized to other modes of ARPACK

// Mult       Functor that computes the matrix-vector multiply in-place, given an array of complex<double>
// n          vector size
// NumEigen   number of eigenvalues to calculate
// tol        eigensolver tolerance
// OutVectors pointer to a vector that stores the eigenvectors.  If you don't need eigenvectors, set to nullptr
// ncv        length of the Krylov sequence.  Must be > NumEigen, recommend at least NumEigen*2.  If zero,
//            then initialize to a default (2*NumEigen+10)
// Sort       if true, then sort the eigenvalues in order of decreasing magnitude
// Verbose    verbose output level

enum class WhichEigenvalues { LargestMagnitude, SmallestMagnitude, LargestReal, SmallestReal,
      LargestImag, SmallestImag, LargestAlgebraic, SmallestAlgebraic, BothEnds };

std::string ToStr(WhichEigenvalues);

// For generic problems, WhichEigenvalues must be one of
// LargestMagnitude, SmallestMagnitude, LargestReal, SmallestReal, LargestImag, SmallestImag
// if InitialGuess is nullptr, then let ARPACK choose the initial residual randomly
template <typename MultFunc>
Vector<std::complex<double>>
DiagonalizeARPACK(MultFunc Mult, int n, int NumEigen, WhichEigenvalues which, std::complex<double> const* InitialGuess, double tol = 1e-10, std::vector<std::complex<double>>* OutputVectors = NULL, int ncv = 0, bool Sort = false, int Verbose = 0);

// Get the eigenvalues with the largest magnitude
template <typename MultFunc>
Vector<std::complex<double>>
DiagonalizeARPACK(MultFunc Mult, int n, int NumEigen, std::complex<double> const* InitialGuess, double tol = 1e-10,
                  std::vector<std::complex<double>>* OutputVectors = NULL,
                  int ncv = 0, bool Sort = false, int Verbose = 0)
{
   return DiagonalizeARPACK(Mult, n, NumEigen, WhichEigenvalues::LargestMagnitude, InitialGuess, tol,
                  OutputVectors, ncv, Sort, Verbose);
}

} // namespace LinearAlgebra

#include "arpack_wrapper.cc"

#endif
