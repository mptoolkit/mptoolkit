// -*- C++ -*- $Id$

#include "linearalgebra/matrix.h"
#include <iostream>
#include <iomanip>

using namespace LinearAlgebra;
using namespace std;

int main()
{
  Matrix<std::complex<double> > M(3,4,std::complex<double>(2.0,1.0));
  Matrix<std::complex<double> > N(3,4,4.0);
  Matrix<std::complex<double> > R; R = M*transpose(conj(N));
  std::cout << R;
}
