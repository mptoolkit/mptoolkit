// -*- C++ -*- $Id$

#include "linearalgebra/matrix.h"
#include "linearalgebra/eigen.h"
#include "linearalgebra/matrix_utility.h"
#include "linearalgebra/fixedvector.h"

#include <iostream>
#include <iomanip>

using namespace LinearAlgebra;

typedef std::complex<double> complex;

template <typename Mat>
void ZeroUpper(Mat& M)
{
   for (unsigned i = 0; i < M.size1(); ++i)
   {
      for (unsigned j = i+1; j < M.size2(); ++j)
      {
	 M(i,j) = 0.0;
      }
   }
}

template <typename Mat>
void ZeroLower(Mat& M)
{
   for (unsigned j = 0; j < M.size2(); ++j)
   {
      for (unsigned i = j+1; i < M.size1(); ++i)
      {
	 M(i,j) = 0.0;
      }
   }
}

int main()
{
   int const dim = 8;

   Matrix<std::complex<double> > I = diagonal_matrix(FixedVector<double>(dim, 1.0));
   Matrix<std::complex<double> > H = random_matrix<std::complex<double> >(dim, dim);
   H = H * herm(H) + I;

   Matrix<std::complex<double> > HChol(H);
   CholeskyFactorizeUpper(HChol);
   ZeroLower(HChol);
   CHECK_CLOSE(herm(HChol) * HChol, H)(HChol);

   HChol = H;
   CholeskyFactorizeLower(HChol);
   ZeroUpper(HChol);
   CHECK_CLOSE(HChol * herm(HChol), H)(HChol);

   Matrix<complex, ColMajor> H2(H);
   CholeskyFactorizeUpper(H2);
   ZeroLower(H2);
   CHECK_CLOSE(herm(H2) * H2, H)(H2);

   H2 = H;
   CholeskyFactorizeLower(H2);
   ZeroUpper(H2);
   CHECK_CLOSE(H2 * herm(H2), H)(H2);
}
