// -*- C++ -*- $Id$

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

template <typename MultFunc>
Vector<std::complex<double> > 
DiagonalizeARPACK(MultFunc Mult, int n, int NumEigen, double tol = 1e-10,
                  std::vector<std::complex<double> >* OutputVectors = NULL,
                  int ncv = 0, bool Sort = false, int Verbose = 0);

// Match left and right eigenvectors to the correct complex eigenvalues.
// This works by finding a pairing of numerically identical eigenvalues, 
// and re-ordering the output arrays so that the corresponding eigenvectors
// are in the same array index.
// If the eigenvalues differ in absolute magnitude by more than tol, then print a warning message to cerr.
// TODO: this doesn't handle degeneracies properly.

void
MatchEigenvectors(int n, 
                  Vector<std::complex<double> >& LeftValues, 
                  std::vector<std::complex<double> >& LeftVectors,
                  Vector<std::complex<double> >& RightValues,
                  std::vector<std::complex<double> >& RightVextors, double tol = 1e-10);

} // namespace LinearAlgebra

#include "arpack_wrapper.cc"

#endif
