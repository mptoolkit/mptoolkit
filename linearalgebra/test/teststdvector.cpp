// -*- C++ -*- $Id$

#include "linearalgebra/stdvector.h"
#include "testdensevector.h"
#include "linearalgebra/vectoroperations.h"
#include "linearalgebra/vectoraddition.h"
#include "linearalgebra/vectortransform.h"
#include "linearalgebra/index.h"
#include "linearalgebra/scalar.h"
#include "linearalgebra/vector.h"

using namespace LinearAlgebra;


void subset(Vector<int> const& v, Vector<int> const& selection)
 {
 Vector<int> res;
 noalias(res)  = v[selection];
 }


int main()
{
   double y = LinearAlgebra::norm_1(-2.1);
   y = LinearAlgebra::norm_2(-2.1);

   test_dense_vector<std::vector<double> >();
   test_real_vector<std::vector<double> >();

   std::vector<double> Vec(10);
   assign(Vec, LinearAlgebra::Range(10,20));
   test_minmax(Vec);
}
