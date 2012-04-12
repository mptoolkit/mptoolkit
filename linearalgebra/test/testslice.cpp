// -*- C++ -*- $Id$

#include "linearalgebra/vector.h"
#include "linearalgebra/slice.h"

int main()
{
   LinearAlgebra::Slice s(0,5,4);
   LinearAlgebra::Vector<double> v1(s);
   std::cout << s << std::endl;
   std::cout << v1 << std::endl;
}
