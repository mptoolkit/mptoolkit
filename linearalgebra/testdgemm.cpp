// -*- C++ -*- $Id$

#include "densematrix.h"

int main()
{
   DenseMatrix<double> Result(4,5,0.0);
   DenseMatrix<double> M1(4,6,1.0);
   DenseMatrix<double> M2(6,5,2.0);
   for (int i = 0; i < M1.rows(); ++i)
   {
      for (int j = 0; j < M1.cols(); ++j)
      {
	 M1(i,j) = i-j;
      }
   }

   for (int i = 0; i < M2.rows(); ++i)
   {
      for (int j = 0; j < M2.cols(); ++j)
      {
	 M2(i,j) = i*j;
      }
   }


   add_product_AB(Result, M1, M2);

   BLAS::DGS::DebugPrintDgemm();
}
