// -*- C++ -*- $Id$

#include "linearalgebra/matrix.h"

using namespace LinearAlgebra;

int main()
{
   LinearAlgebra::Matrix<double> M(3,3, 2.0);
   LinearAlgebra::Matrix<double> N(3,3, 3.0);
   
   LinearAlgebra::Matrix<double> R;

   R = transform(M, N, LinearAlgebra::Multiplication<double, double>());

   TRACE(R);

   LinearAlgebra::Matrix<double, RowMajor> MM(2,3, 2.0);
   LinearAlgebra::Matrix<double, ColMajor> NN(2,3, 3.0);

   TRACE(transform(MM, NN, LinearAlgebra::Multiplication<double, double>()));

   TRACE(MM+NN);

   TRACE(element_prod(MM,NN));
}

