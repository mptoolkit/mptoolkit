// -*- C++ -*- $Id$

#include "linearalgebra/hashvector.h"
#include "testvectorsparse.h"
#include "linearalgebra/matrix.h"
#include "linearalgebra/matrix.h"
#include "linearalgebra/matrixmatrixmultiplication.h"

int main()
{
   test_real_sparse_vector<LinearAlgebra::HashVector<double> >();

#if 0
   LinearAlgebra::HashVector<LinearAlgebra::Matrix<double> > v1(10), v2(10);
   TRACE(typeid(parallel_prod(v1,v2)).name());
#endif
}
