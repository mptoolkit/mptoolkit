// -*- C++ -*- $Id$

#include "linearalgebra/matrix.h"
#include "linearalgebra/vector.h"

using namespace LinearAlgebra;

int main()
{
   Matrix<double> M(10,10);

   fill(M, 2.0);
   CHECK_CLOSE(norm_frob(M), 20);

   Matrix<std::complex<double> > X(10, 10);
   fill(real(X), 3.0);
   fill(imag(X), 4.0);
   CHECK_CLOSE(norm_frob(X), 50);

   Vector<double> V(100);
   fill(V, 2.0);
   CHECK_CLOSE(norm_frob(V), 20);

   Vector<std::complex<double> > W(100);
   fill(real(W), 3.0);
   fill(imag(W), 4.0);
   CHECK_CLOSE(norm_frob(W), 50);
}
