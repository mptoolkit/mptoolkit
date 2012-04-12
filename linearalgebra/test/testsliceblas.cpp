// -*- C++ -*- $Id$


#define BLAS1_TRACE_DETAILED

#include "linearalgebra/matrix.h"
#include "linearalgebra/vector.h"
#include "common/trace.h"

using namespace LinearAlgebra;

using tracer::typeid_name;

int main()
{
   Vector<double> V(100);
   Vector<double> X(100);

   X[Slice(1,1,1)] = V[Slice(1,1,1)];

   Matrix<double> M(100,100);
   matrix_row(M, 10)[Slice(1,1,1)] = V[Slice(1,1,1)];
}
