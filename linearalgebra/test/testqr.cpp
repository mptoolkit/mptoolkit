// -*- C++ -*- $Id$

#include "linearalgebra/matrix.h"
#include "linearalgebra/eigen.h"
#include "linearalgebra/matrix_utility.h"
#include "linearalgebra/fixedvector.h"

#include <iostream>
#include <iomanip>

using namespace LinearAlgebra;

typedef std::complex<double> complex;

void Test(int m, int n)
{
   Matrix<std::complex<double> > H = random_matrix<std::complex<double> >(m, n);
   Matrix<std::complex<double> > R(H);
   Matrix<std::complex<double> > Q = QR_Factorize(R);
   CHECK_CLOSE(H,Q*R)(Q)(R)(conj(Q)*R)(herm(Q)*R)(trans(Q)*R);
}

int main()
{
   Test(2,2);
   Test(3,3);
   Test(4,4);
   Test(2,1);
   Test(2,3);
   Test(10,5);
   Test(5,10);
}
