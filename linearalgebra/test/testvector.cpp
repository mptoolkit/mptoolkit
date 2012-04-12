// -*- C++ -*- $Id$

#include "linearalgebra/vector.h"
#include "testdensevector.h"

int main()
{
   test_dense_vector<LinearAlgebra::Vector<double> >();

   test_real_vector<LinearAlgebra::Vector<double> >();

   typedef LinearAlgebra::Vector<LinearAlgebra::Vector<double> > vt;

   test_ctor<vt>();
   
   LinearAlgebra::Vector<double> v1(4,1);
   test_dense_single<vt>(v1);

   vt v2(10, v1);
   test_assign(v2);

   test_double_negate<vt>();
   test_real_scalar_vector_nop<vt>(v2);
   test_equal_to(v2);

   LinearAlgebra::Vector<double> V(1);
   V = V[LinearAlgebra::Range(0,0)];
   CHECK_EQUAL(size(V), 0);
}
