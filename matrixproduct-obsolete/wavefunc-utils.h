// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// matrixproduct-obsolete/wavefunc-utils.h
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

#if !defined(WAVEFUNC_UTILS_74YT743YT7YEIDSU)
#define WAVEFUNC_UTILS_74YT743YT7YEIDSU

#include "matrixproduct/centerwavefunction.h"
#include "matrixproduct/splitoperator.h"
#include "matrixproduct/operatorstack.h"

typedef OperatorStack<MPStateComponent> SuperblockOperator;
typedef OperatorStack<MatrixOperator>   TransformOperator;

// given two wavefunctions x,y, at the same rotation and with the same quantum number,
// initializes the TransformOperator that maps between x and y.
void InitializeTransformStack(TransformOperator& x_y, 
                              CenterWavefunction const& x, 
                              CenterWavefunction const& y);

// Update the transform matrix elements, after x and y have been rotated
void TransformStackRotateLeft(TransformOperator& x_y, 
			      CenterWavefunction const& x, 
			      CenterWavefunction const& y);

void TransformStackRotateRight(TransformOperator& x_y, 
			       CenterWavefunction const& x, 
			       CenterWavefunction const& y);

// Update the transform matrix elements, after x and y have been modified
void TransformStackUpdateLeft(TransformOperator& x_y, 
			      CenterWavefunction const& x, 
			      CenterWavefunction const& y);

void TransformStackUpdateRight(TransformOperator& x_y, 
			       CenterWavefunction const& x, 
			       CenterWavefunction const& y);

// Initializes the matrix elements x_A_y.
// Also rotates A to be the correct location, as necessary.
void InitializeSuperblockStack(SuperblockOperator& x_A_y, 
			       CenterWavefunction const& x, 
			       SplitOperator& A,
			       CenterWavefunction const& y);

// Update the transform matrix elements, after x and y have been rotated
void SuperblockStackRotateLeft(SuperblockOperator& x_A_y, 
			       CenterWavefunction const& x, 
			       SplitOperator const& A,
			       CenterWavefunction const& y);
void SuperblockStackRotateRight(SuperblockOperator& x_A_y, 
				CenterWavefunction const& x, 
				SplitOperator const& A,
				CenterWavefunction const& y);

// Update the transform matrix elements, after x and y have been updated
void SuperblockStackUpdateLeft(SuperblockOperator& x_A_y, 
			       CenterWavefunction const& x, 
			       SplitOperator const& A,
			       CenterWavefunction const& y);
void SuperblockStackUpdateRight(SuperblockOperator& x_A_y, 
				CenterWavefunction const& x, 
				SplitOperator const& A,
				CenterWavefunction const& y);

#endif
