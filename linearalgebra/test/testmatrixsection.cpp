// -*- C++ -*- $Id$

#include "linearalgebra/matrix.h"
#include "linearalgebra/matrixsection.h"
#include "linearalgebra/matrixtranspose.h"

using namespace LinearAlgebra;

int main()
{
   Matrix<double> M(3,3,0);
   M(0,0) = 1;
   M(0,1) = 5;
   M(0,2) = 4;
   M(1,0) = 3;
   M(1,1) = 8;
   M(1,2) = 7;
   M(2,0) = 9;
   M(2,1) = 6;
   M(2,2) = 2;

   TRACE(M);

   Range r(0,2);
   Slice s(0,2,2);


   BOOST_MPL_ASSERT((is_proxy_reference<Range const&>));
   BOOST_MPL_ASSERT((is_const_proxy_reference<Range const&>));

   MatrixSection<Matrix<double>&, Range const&, Slice const&> MySec(M, r, s);
   TRACE(MySec);

   get_element(MySec, 1,1) = -10;
   TRACE(MySec);
   TRACE(M);

   Matrix<double> Test(MySec);

   MySec = transpose(Test);

   TRACE(Test)(MySec)(M);

   TRACE(project(M, r, s));

   TRACE(M);
   M(r,s) *= 2.0;
   TRACE(M);
}
