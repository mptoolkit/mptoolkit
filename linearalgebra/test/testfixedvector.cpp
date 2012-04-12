// -*- C++ -*- $Id$

#include "linearalgebra/fixedvector.h"
#include "linearalgebra/scalar.h"
#include <cmath>

using namespace LinearAlgebra;

int main()
{
   FixedVector<double> v(10, 2.0);
   CHECK_EQUAL(norm_1(v), 20.0);

   FixedVector<double> c(v);
   CHECK_EQUAL(c,v);

   multiply(v, 2);
   CHECK_EQUAL(norm_1(v), 40.0);

   FixedVector<std::complex<double> > w(10, std::complex<double>(0,4.0));
   CHECK_EQUAL(norm_1(v), 40.0);

   FixedVector<std::complex<double> > y = v+w;
   CHECK_CLOSE(norm_1(y), 10 * std::sqrt(32.0));

   FixedVector<std::complex<double> > z;
   z = v+w*2.0;
   CHECK_CLOSE(norm_1(z), 10 * std::sqrt(80.0));

   real(z) = FixedVector<double>(10, 1.0);
   CHECK_CLOSE(norm_1(z), 10 * sqrt(65.0));
}
