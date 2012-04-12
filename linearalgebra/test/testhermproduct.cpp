// -*- C++ -*- $Id$

#include "linearalgebra/matrix.h"

using namespace LinearAlgebra;

typedef std::complex<double> complex;

int main()
{
   Matrix<complex> M1(2,3,1.0);
   Matrix<complex> M2(3,4,2.0);
   
   TRACE(herm(M1*M2));
}
