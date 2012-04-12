// -*- C++ -*- $Id$

#include "linearalgebra/matrix.h"

using namespace LinearAlgebra;

int main()
{
   Matrix<double> M(4,4,0.0);
   TRACE( M(range(0,3),all)(range(0,2),all) );

   TRACE(  M(range(0,3),all)(range(0,2),all)(1,1) );
}
