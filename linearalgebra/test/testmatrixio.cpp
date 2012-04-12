// -*- C++ -*- $Id$

#include "linearalgebra/matrix.h"
#include "linearalgebra/sparsematrix.h"

using namespace LinearAlgebra;

int main()
{
   Matrix<double> M(10,10,0.1);

   std::cout << M(all,range(1,5))(5,all) << '\n';
}
