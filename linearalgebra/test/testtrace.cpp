// -*- C++ -*- $Id$

#include "linearalgebra/matrix.h"

using namespace LinearAlgebra;

int main()
{
  Matrix<double> A(2,2);
  A(0,0) = 3;
  A(1,1) = -1;
  A(1,0) = 2;
  A(0,1) = 8;
  CHECK_EQUAL(trace(A), 2);
  CHECK_EQUAL(inner_prod(A,A), 78);
  CHECK_EQUAL(inner_prod(A,trans(A)), 42);
}
