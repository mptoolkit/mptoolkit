// -*- C++ -*- $Id$

#include "linearalgebra/matrix.h"
#include "linearalgebra/eigen.h"
#include "linearalgebra/matrix_utility.h"

#include <iostream>
#include <iomanip>

using namespace LinearAlgebra;

int main()
{
   int const dim = 8;
   Matrix<double> H = random_matrix<double>(dim, dim);
   H=H+transpose(H);

   // simplest case
   Matrix<double> HCopy(H);
   Vector<double> Eigen = DiagonalizeSymmetric(H);
   CHECK_CLOSE(H * HCopy * transpose(H), diagonal_matrix(Eigen));

   Matrix<double> T = H;
   H = HCopy;

   // column major matrix
   Matrix<double, ColMajor> HCol = H;
   Vector<double> Eigen2 = DiagonalizeSymmetric(HCol);
   CHECK_CLOSE(Eigen2, Eigen);
   CHECK_CLOSE(HCol, T);

   // a slice expression
   Vector<double> EigenSlice = DiagonalizeSymmetric(H(Slice(2,3,2), Slice(2,3,2)));
   Matrix<double> HSlice = H(Slice(2,3,2), Slice(2,3,2));
   Matrix<double> HCopySlice = HCopy(Slice(2,3,2), Slice(2,3,2));
   CHECK_CLOSE(HSlice * HCopySlice * transpose(HSlice), diagonal_matrix(EigenSlice));

   // test on a non-trivial expression (real & imaginary parts of a complex matrix)
   Matrix<std::complex<double> > Z = random_matrix<std::complex<double> >(dim, dim);
   Z=Z+transpose(Z);
   Matrix<std::complex<double> > ZCopy = Z;
   Vector<double> ZRealEigen = DiagonalizeSymmetric(real(Z));
   CHECK_CLOSE(real(Z) * real(ZCopy) * transpose(real(Z)), diagonal_matrix(ZRealEigen));

   Vector<double> ZImagEigen = DiagonalizeSymmetric(imag(Z));
   CHECK_CLOSE(imag(Z) * imag(ZCopy) * transpose(imag(Z)), diagonal_matrix(ZImagEigen));
}
