// -*- C++ -*- $Id$

#include "linearalgebra/vector.h"

using namespace LinearAlgebra;

int main()
{
  Vector<double> v(3);
  v[0] = -10;
  v[1] = 3;
  v[2] = 7;
  CHECK_EQUAL(norm_inf(v), 10);

  Vector<std::complex<double> > vc(10, std::complex<double>(3,4));
  CHECK_EQUAL(norm_inf(vc), 5);
}
