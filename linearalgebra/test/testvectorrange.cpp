// -*- C++ -*- $Id$

#include "linearalgebra/vector.h"

using namespace LinearAlgebra;

int main()
{
   Vector<int> v = range(0,5);
   TRACE(v);
   v = v[range(2,5)];
   TRACE(v);
}
