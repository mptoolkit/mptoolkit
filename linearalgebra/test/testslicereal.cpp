// -*- C++ -*- $Id$

#include "linearalgebra/matrix.h"
#include <complex>

using namespace LinearAlgebra;

int main()
{
   Matrix<std::complex<double> > M(10,10);

   real(M(all, 1)) *= 2;
   imag(M(all, 1)) *= 2;
}
