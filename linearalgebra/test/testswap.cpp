// -*- C++ -*- $Id$

#include "linearalgebra/vector.h"
#include "linearalgebra/matrix.h"

using namespace LinearAlgebra;

int main()
{
   Vector<double> V(2);
   V[0] = 3;
   V[1] = 4;

   Vector<double> V2(3);
   V2[0] = -3;
   V2[1] = -4;
   V2[2] = 5;

   Vector<double> U = V;
   Vector<double> U2 = V2;

   swap(V, V2);
   CHECK_EQUAL(V, U2);
   CHECK_EQUAL(V2, U);
}
