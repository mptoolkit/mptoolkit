// -*- C++ -*- $Id$

#include "linearalgebra/matrix.h"
#include <iostream>
#include <iomanip>

using namespace LinearAlgebra;
using namespace std;

int main()
{
  Matrix<double> M(3,4,2.0);
  Matrix<double> N(3,4,4.0);
  M(0,1) = 1;

  Matrix<double> MNCheck(3,4,8.0);
  MNCheck(0,1) = 4;

  Matrix<double> MNElementProd = apply_func(BinaryOperator<Multiplication, double, double>(), M, N);

  CHECK_EQUAL(MNElementProd, MNCheck);
}
