// -*- C++ -*- $Id$

#include "linearalgebra/matrix.h"
#include "linearalgebra/sparsematrix.h"
#include "linearalgebra/matrix_utility.h"
#include <complex>

using namespace LinearAlgebra;

int main()
{
   typedef Matrix<std::complex<double> > T;
   T M1 = random_matrix<T::value_type>(2,2);
   SparseMatrix<T> M(4,5);
   M(0,2) = M1;
   TRACE(conj(-M));
   TRACE(herm(M));
   TRACE(transform(M, Herm<T>()));
   SparseMatrix<T> N(herm(M));
}
