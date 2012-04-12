// -*- C++ -*- $Id$

#include "linearalgebra/vector.h"
#include "linearalgebra/matrix.h"
#include "linearalgebra/scalar.h"

using namespace LinearAlgebra;

int main()
{
   Vector<double> v = Range(1,10);
   CHECK_EQUAL(sum(v), 45.0);

   Matrix<double> M(3,3,1.0);
   Vector<Matrix<double> > vm(3, M);

   Matrix<double> N(3,3,0.0);

   typedef boost::mpl::print<Sum<Vector<Matrix<double> > >::result_type>::type d;


   N += sum(vm);
   CHECK_EQUAL(N, 3*M);

   Vector<Matrix<double> > vn;
   N += sum(vn);
   CHECK_EQUAL(N, 3*M);
}
