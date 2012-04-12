// -*- C++ -*- $Id$

#include "linearalgebra/matrix.h"

using namespace LinearAlgebra;

int main()
{
   Matrix<double> M(2,3);
   M(0,0) = 1;   M(0,1) = 2;   M(0,2) = 3;
   M(1,0) = 4;   M(1,1) = 5;   M(1,2) = 6;

   Matrix<double> V(1,2);
   V(0,0) = -1;
   V(0,1) = 1;

   Matrix<double> V2 = transpose(V);

   Matrix<double> W(3,1);
   W(0,0) = -2;
   W(1,0) = 3;
   W(2,0) = -4;

   Matrix<double> W2 = transpose(W);

   Matrix<double, RowMajor> TempR;
   Matrix<double, ColMajor> TempC;

   Matrix<double> MWCheck(2,1);
   MWCheck(0,0) = -8;
   MWCheck(1,0) = -17;

   TempR = M*W;
   CHECK_CLOSE(MWCheck, TempR);
   TempR = transpose(W) * transpose(M);
   CHECK_CLOSE(transpose(MWCheck), TempR);
   TempR = W2 * transpose(M);
   CHECK_CLOSE(transpose(MWCheck), TempR);

   TempC = M*W;
   CHECK_CLOSE(MWCheck, TempC);
   TempC = transpose(W) * transpose(M);
   CHECK_CLOSE(transpose(MWCheck), TempC);
   TempC = W2 * transpose(M);
   CHECK_CLOSE(transpose(MWCheck), TempC);

   TempR = -MWCheck;
   TempR += M*W;
   CHECK_CLOSE(norm_frob(TempR), 0);
   TempR = -transpose(MWCheck);
   TempR += transpose(W) * transpose(M);
   CHECK_CLOSE(norm_frob(TempR), 0);
   TempR = -transpose(MWCheck);
   TempR += W2 * transpose(M);
   CHECK_CLOSE(norm_frob(TempR), 0);

   TempC = -MWCheck;
   TempC += M*W;
   CHECK_CLOSE(norm_frob(TempC), 0);
   TempC = -transpose(MWCheck);
   TempC += transpose(W) * transpose(M);
   CHECK_CLOSE(norm_frob(TempC), 0);
   TempC = -transpose(MWCheck);
   TempC += W2 * transpose(M);
   CHECK_CLOSE(norm_frob(TempC), 0);

   TempR = MWCheck;
   TempR -= M*W;
   CHECK_CLOSE(norm_frob(TempR), 0);
   TempR = transpose(MWCheck);
   TempR -= transpose(W) * transpose(M);
   CHECK_CLOSE(norm_frob(TempR), 0);
   TempR = transpose(MWCheck);
   TempR -= W2 * transpose(M);
   CHECK_CLOSE(norm_frob(TempR), 0);

   TempC = MWCheck;
   TempC -= M*W;
   CHECK_CLOSE(norm_frob(TempC), 0);
   TempC = transpose(MWCheck);
   TempC -= transpose(W) * transpose(M);
   CHECK_CLOSE(norm_frob(TempC), 0);
   TempC = transpose(MWCheck);
   TempC -= W2 * transpose(M);
   CHECK_CLOSE(norm_frob(TempC), 0);

   Matrix<double> VMCheck(1,3);
   VMCheck(0,0) = 3;
   VMCheck(0,1) = 3;
   VMCheck(0,2) = 3;

   TempR = V*M;
   CHECK_CLOSE(VMCheck, TempR);
   TempR = transpose(M) * transpose(V);
   CHECK_CLOSE(transpose(VMCheck), TempR);
   TempR = transpose(M) * V2;
   CHECK_CLOSE(transpose(VMCheck), TempR);

   TempC = V*M;
   CHECK_CLOSE(VMCheck, TempC);
   TempC = transpose(M) * transpose(V);
   CHECK_CLOSE(transpose(VMCheck), TempC);
   TempC = transpose(M) * V2;
   CHECK_CLOSE(transpose(VMCheck), TempC);

   TempR = -VMCheck;
   TempR += V*M;
   CHECK_CLOSE(norm_frob(TempR), 0);
   TempR = -transpose(VMCheck);
   TempR += transpose(M) * transpose(V);
   CHECK_CLOSE(norm_frob(TempR), 0);
   TempR = -transpose(VMCheck);
   TempR += transpose(M) * V2;
   CHECK_CLOSE(norm_frob(TempR), 0);

   TempC = -VMCheck;
   TempC += V*M;
   CHECK_CLOSE(norm_frob(TempC), 0);
   TempC = -transpose(VMCheck);
   TempC += transpose(M) * transpose(V);
   CHECK_CLOSE(norm_frob(TempC), 0);
   TempC = -transpose(VMCheck);
   TempC += transpose(M) * V2;
   CHECK_CLOSE(norm_frob(TempC), 0);

   TempR = VMCheck;
   TempR -= V*M;
   CHECK_CLOSE(norm_frob(TempR), 0);
   TempR = transpose(VMCheck);
   TempR -= transpose(M) * transpose(V);
   CHECK_CLOSE(norm_frob(TempR), 0);
   TempR = transpose(VMCheck);
   TempR -= transpose(M) * V2;
   CHECK_CLOSE(norm_frob(TempR), 0);

   TempC = VMCheck;
   TempC -= V*M;
   CHECK_CLOSE(norm_frob(TempC), 0);
   TempC = transpose(VMCheck);
   TempC -= transpose(M) * transpose(V);
   CHECK_CLOSE(norm_frob(TempC), 0);
   TempC = transpose(VMCheck);
   TempC -= transpose(M) * V2;
   CHECK_CLOSE(norm_frob(TempC), 0);
}
