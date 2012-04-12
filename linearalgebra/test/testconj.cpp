// -*- C++ -*- $Id$

#include "linearalgebra/matrix.h"
#include "linearalgebra/vector.h"
#include <iostream>
#include <iomanip>

using namespace LinearAlgebra;
using namespace std;

int main()
{
  double y = 3;
  double x = conj(y);
  CHECK_EQUAL(x,y);

  Vector<double> MV(3,-1);
  Vector<double> NV = conj(MV);
  CHECK_CLOSE(MV, NV);

  Vector<complex<double> > ZV(3,complex<double>(0,1));
  Vector<complex<double> > CV(conj(ZV));
  CHECK_CLOSE(CV, -ZV);

  Matrix<double> M(3,3,-1);
  Matrix<double> N = conj(M);
  CHECK_CLOSE(M, N);

  Matrix<complex<double> > Z(3,3,complex<double>(0,1));
  Matrix<complex<double> > C(conj(Z));
  CHECK_CLOSE(C, -Z);
}
