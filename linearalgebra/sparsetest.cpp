// -*- C++ -*- $Id$

#include "sparsematrix.h"

using namespace LinearAlgebra;

int main()
{
  SparseMatrix<double> M(20, 20);

  M(0,0) = 30;
  M(0,1) = 40;
  M(1,0) = -40;

  std::cout << M << std::endl;
}

