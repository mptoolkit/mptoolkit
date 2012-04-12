// -*- C++ -*- $Id$

#include "linearalgebra/matrix.h"

using namespace LinearAlgebra;

int main()
{
   Matrix<double> A, B;
   Matrix<double> M2(2,2,1);
   Matrix<double> M3(3,3,1);
   Matrix<double> N(2,2,1E-10);
   Matrix<double> Z(2,2,0);

   CHECK(equal(A, B));
   CHECK(equal(M2, M2));
   CHECK(equal(M2+M2, M2*2));
   CHECK(equal(M2,M2+Z));
   CHECK(!equal(M2,M3));
   CHECK(!equal(M2,M2+N));
}
