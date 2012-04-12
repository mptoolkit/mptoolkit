// -*- C++ -*- $Id$

#include "linearalgebra/matrix.h"
#include "linearalgebra/matrix_utility.h"

using namespace LinearAlgebra;

int main()
{
   Matrix<double> M = random_matrix<double>(10, 12);

   M *= random_matrix<double>(12, 12);
   M *= random_matrix<double>(12,8);
}
