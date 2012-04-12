// -*- C++ -*- $Id$

#include "linearalgebra/vector.h"
#include <string>

using namespace LinearAlgebra;

int main()
{
   Vector<std::string> X(3);
   X[0] = "hello";
   X[1] = "world";
   X[2] = "goodbye";

   X = X[range(1,2)];
   CHECK_EQUAL(size(X), 1);
   CHECK_EQUAL(X[0], "world");
}
