// -*- C++ -*- $Id$

#include "linearalgebra/matrixdirectproduct.h"
#include "linearalgebra/matrixsection.h"
#include "linearalgebra/matrixtranspose.h"
#include "linearalgebra/matrixmatrixmultiplication.h"
#include "linearalgebra/matrix.h"
#include <iostream>
#include <iomanip>

using namespace LinearAlgebra;

int main()
{
   Matrix<double> M1(2,2,0.0), M2(3,3,0.0);

   M1(0,0) = 1;
   M1(0,1) = 2;
   M1(1,0) = 3;
   M1(1,1) = 4;

   M2(0,0) = -1;
   M2(0,1) = -5;
   M2(0,2) = 3;
   M2(1,0) = 7;
   M2(1,1) = 2;
   M2(1,2) = 8;
   M2(2,0) = 0;
   M2(2,1) = 2;
   M2(2,2) = 1;

   TRACE(M1)(M2)(direct_product(M1,M2));

   Matrix<Matrix<double> > MM1(1,1,M1), MM2(1,1,M2);

   TRACE(direct_product(MM1, M2));

   TRACE(direct_product(M1, MM2));

   Matrix<Matrix<double> > MM(direct_product(MM1, MM2));

   TRACE(MM);

   TRACE(direct_product(MM1, MM2));

}
