// -*- C++ -*- $Id$

#include "linearalgebra/matrix.h"
#include "linearalgebra/vector.h"

using namespace LinearAlgebra;

int main()
{
   Matrix<double> M(2,2);
   M(0,0) = 5;
   M(0,1) = 10;
   M(1,0) = -3;
   M(1,1) = 7;

   Vector<double> FirstRow(2);
   FirstRow[0] = 5;
   FirstRow[1] = 10;

   Vector<double> SecondRow(2);
   SecondRow[0] = -3;
   SecondRow[1] = 7;

   Vector<double> FirstCol(2);
   FirstCol[0] = 5;
   FirstCol[1] = -3;

   Vector<double> SecondCol(2);
   SecondCol[0] = 10;
   SecondCol[1] = 7;

   CHECK_CLOSE(matrix_row(M, 0), FirstRow);
   CHECK_CLOSE(matrix_row(M, 1), SecondRow);

   CHECK_CLOSE(project1(M, 0), FirstRow);
   CHECK_CLOSE(project1(M, 1), SecondRow);

   CHECK_CLOSE(M(0, all), FirstRow);
   CHECK_CLOSE(M(1, all), SecondRow);


   CHECK_CLOSE(matrix_col(M, 0), FirstCol);
   CHECK_CLOSE(matrix_col(M, 1), SecondCol);

   CHECK_CLOSE(project2(M, 0), FirstCol);
   CHECK_CLOSE(project2(M, 1), SecondCol);

   CHECK_CLOSE(M(all, 0), FirstCol);
   CHECK_CLOSE(M(all, 1), SecondCol);


   matrix_row(M, 1) *= 2;

   CHECK_CLOSE(project1(M, 0), FirstRow);
   CHECK_CLOSE(project1(M, 1), 2 * SecondRow);

   project1(M, 1) *= 2;

   CHECK_CLOSE(project1(M, 0), FirstRow);
   CHECK_CLOSE(project1(M, 1), 4 * SecondRow);

   M(1, all) *= 2;

   CHECK_CLOSE(project1(M, 0), FirstRow);
   CHECK_CLOSE(project1(M, 1), 8 * SecondRow);

   M(1, all) = M(0, all) + M(0, all);
   CHECK_CLOSE(M(1, all), 2 * M(0, all));

   M(0,all) += M(1,1)*(M(0,all)+M(1,all));
   CHECK_EQUAL(M(0,0), 305);
   CHECK_EQUAL(M(0,1), 610);
   CHECK_EQUAL(M(1,0), 10);
   CHECK_EQUAL(M(1,1), 20);
}
