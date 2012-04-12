// -*- C++ -*- $Id$

#include "linearalgebra/matrix.h"
#include "linearalgebra/eigen.h"

#include <iostream>
#include <iomanip>

using namespace std;
using namespace LinearAlgebra;

int main()
{
  Matrix<double> A(2,2, 0);
  Matrix<double> B(2,2, 0);
  B(0,0) = 4;  B(1,1) = 4;  B(1,0) = B(0,1) = 1;

  A(0,0) = 1; A(1,1) = 1;
  A(1,0) = A(0,1) = 4;

  Vector<double> Evalues;
  Matrix<double> Evectors;

  GeneralizedEigenHermitian(A, B, Evalues, Evectors, Range(0,2));

  Matrix<double> ExpectedEvectors(2,2);
  ExpectedEvectors(0,0) = -1.0 / std::sqrt(6.0);
  ExpectedEvectors(0,1) = 1.0 / std::sqrt(6.0);
  ExpectedEvectors(1,0) = ExpectedEvectors(1,1) = 1.0 / std::sqrt(10.0);
  Vector<double> ExpectedEvalues(2);
  ExpectedEvalues[0] = -1;
  ExpectedEvalues[1] = 1;

  CHECK(equal(Evalues, ExpectedEvalues));
  CHECK(equal(Evectors, ExpectedEvectors));

}
