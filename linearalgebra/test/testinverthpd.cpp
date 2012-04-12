// -*- C++ -*- $Id$

#include "linearalgebra/matrix.h"
#include "linearalgebra/eigen.h"
#include "linearalgebra/matrix_utility.h"
#include "linearalgebra/fixedvector.h"

#include <iostream>
#include <iomanip>

using namespace LinearAlgebra;

int main()
{
   int const dim = 8;

   Matrix<std::complex<double> > I = diagonal_matrix(FixedVector<double>(dim, 1.0));
   Matrix<std::complex<double> > H = random_matrix<std::complex<double> >(dim, dim);
   H = H * herm(H) + I;

   Matrix<std::complex<double> > HInv(H);

   InvertHPD(HInv);

   CHECK_CLOSE(HInv * H, I)(H)(HInv);
   CHECK_CLOSE(H * HInv, I)(H)(HInv);
}
