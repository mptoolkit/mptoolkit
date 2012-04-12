// -*- C++ -*- $Id$

#include "linearalgebra/eigen.h"
#include "linearalgebra/matrix.h"
#include "linearalgebra/vector.h"
#include "linearalgebra/matrixtranspose.h"
#include "linearalgebra/matrixmatrixmultiplication.h"
#include "linearalgebra/matrixaddition.h"
#include "linearalgebra/scalar.h"
#include "linearalgebra/matrixtransform.h"

using namespace LinearAlgebra;

int main()
{
   Matrix<double> M(10,10);  // left hand side
   Matrix<double, ColMajor> R(10,2);   // 2 right hand side vectors

   for (int i = 0; i < 10; ++i)
   {
      for (int j = 0; j < 10; ++j)
      {
         M(i,j) = i * j + i + j;
         if (i == j) M(i,j) += 1000;
      }

      for (int j = 0; j < 2; ++j)
      {
         R(i,j) = j+i;
      }
   }

   Matrix<double, ColMajor> S = LinearSolveSPD(M,R);
   Matrix<double> RCheck = M * S;

   CHECK(equal(R, RCheck))(R)(RCheck);
}
