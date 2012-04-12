// -*- C++ -*- $Id$

#include "linearalgebra/matrix.h"

using namespace LinearAlgebra;

int main()
{
   Matrix<double> N(4,4,1);
   Matrix<Matrix<double> > NN(4,4,N);
   TRACE( parallel_prod( N, NN ).get() );
   Matrix<double> X = parallel_prod( N, NN );
}
