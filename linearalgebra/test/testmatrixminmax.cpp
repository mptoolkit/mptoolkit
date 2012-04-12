// -*- C++ -*- $Id$

#include "linearalgebra/matrix.h"

using namespace LinearAlgebra;

int main()
{
   Matrix<double> M(3,3,0.0);
   M(0,1) = 1;
   M(1,2) = -5;

   CHECK_EQUAL(min(M), -5);
   CHECK_EQUAL(max(M), 1);

   CHECK_EQUAL(*iter_matrix_max(iterate(M-M)), 0);
}
